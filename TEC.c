#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <omp.h>
#include "utils.h"


#define GF1 1575420000.0
#define GF2 1227600000.0
#define RF1 1602000000.0
#define RF2 1246000000.0
#define DRF1    562500.0
#define DRF2    437500.0
#define H_ION 450000.0
#define SPEED_OF_LIGHT 299792458.0
#define PI  3.1415926535
#define UGLE 0.000072921151467
#define MU 398600441800000.0
#define J20 1082.62982126e-6
#define Ra 6378137.0
#define Rb 6356752.3142
#define RC 6371000.0
#define K_ION 40.308


//структура содержимого n-файла на конкретный момент времени
typedef struct {
	GFile g;
	NFile n;
} navFile; 

typedef struct {
	double X;
	double Y;
	double Z;
	double P;
	double P1;
	double P2;
	double L1;
	double L2;
	double C1;
	double C2;
	double C;
	double dT;
	int K;
	bool satVisible;
} SatParams;//параметры КА необходимы для рассчёта ПЭС

typedef struct {
	Date time;
	char* type;
	SatParams sat[max(GNUM, RNUM)+1];
} Sat;//структура положений КА на конкретный момент времени

VECTOR(SatParams, SatParams)

VECTOR(Sat, Sat)

typedef struct {
	OFileVector o;
	NFileVector n;
	GFileVector g;
} navFiles;//структура содержимого всех файлов

typedef struct {
	NFileVector n;
	GFileVector g;
} prevFile;//структура содержимого n и g файлов на предыдущий день
//вектор строк



//можель движения в относительной неинерциальной системе отсчёта, с учётом переносной, Кориолисовой силы и полярного сжатия Земли
void model(double t, double q[6], double dq[6]) {
	const double r = SQRT(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]);
	
	dq[0] = q[3];
	dq[1] = q[4];
	dq[2] = q[5];
	dq[3] = -MU * q[0] / (r*r*r) + q[0] * UGLE*UGLE + 2 * q[4] * UGLE - 1.5 * J20 * MU * Ra * Ra * (1 - 5 * pow(q[2] / r, 2)) / pow(r, 5)*q[0];
	dq[4] = -MU * q[1] / (r*r*r) + q[1] * UGLE*UGLE - 2 * q[3] * UGLE - 1.5*J20*MU*Ra*Ra*(1 - 5 * pow(q[2] / r, 2)) / pow(r, 5)*q[1];
	dq[5] = -MU * q[2] / (r*r*r) - 1.5*J20 * MU * Ra * Ra * (3 - 5 * pow(q[2] / r, 2)) / pow(r, 5)*q[2];
}

//трилинейная интерполяция карт ПЭС по времени, долготе и широте
double TECmapSpher(IonexFileVector *iF, Date time, double lat,double lon){
	lat *= 180 / PI;
	lon *= 180 / PI;
	if (lat > 87.5 || lat < -87.5) return 9999;// если выходит за пределы карты возвращаем 9999, потом при рассчёте проверяем на это значение и пропускаем если такой результат

	double timeSec = epoch2time(time);
	
	int i,j,k;
	//находим узловые точки по времени
	for ( i = 0; i < iF->size; i++ ) 
		if (timeSec < epoch2time(iF->container[i].time)) break;

	//по широте
	for ( j = 0; j < 71; j++ ) 
		if (lat > iF->container[i].lat[j]) break;
	
	//по долготе
	for ( k = 0; k < 73; k++ ) 
		if (lon < iF->container[i].lon[k]) break;
	
	//интерполируем по 8-ми узловым точкам
	double latKn[2] = { iF->container[i].lat[j-1] ,iF->container[i].lat[j] };
	double lonKn[2] = { iF->container[i].lon[k-1] ,iF->container[i].lon[k] };
	double timeKn[2] = { epoch2time(iF->container[i - 1].time),epoch2time(iF->container[i].time) };
	double TECKn[2][2][2] = { 
		{
			{iF->container[i - 1].TEC[j - 1][k - 1],iF->container[i - 1].TEC[j][k-1]},
			{iF->container[i - 1].TEC[j - 1][k ],iF->container[i - 1].TEC[j ][k]}
		},
		{
			{iF->container[i ].TEC[j - 1][k - 1],iF->container[i ].TEC[j][k - 1]},
			{iF->container[i ].TEC[j - 1][k],iF->container[i].TEC[j][k]}
		}
	};

	for (int _i = 0; _i < 2; _i++) 
		for (int _j = 0; _j < 2; _j++) 
			for (int _k = 0; _k < 2; _k++) 
				if (TECKn[_i][_j][_k] == 9999) return 9999;
			
	return triLinInter(latKn, lonKn, timeKn, TECKn, lat, lon, timeSec) / 10;
}

//рассчёт ПЭС по прямоугольным координатам
double TECmapDec(IonexFileVector *iF, Date time, double X, double Y, double Z) {
	double lat, lon;
	DecToSpher(X, Y, Z, &lat, &lon);
	return TECmapSpher(iF, time, lat, lon);
}


SatParams GSatPosition(NFile ef, int ind, Date date) {
	double n0 = sqrt(MU) / pow(ef.G[ind].A, 3);

	double Tem = time2gpst(epoch2time(date), NULL);
	
	double t = Tem - ef.G[ind].Toe;
	double nn = n0 + ef.G[ind].delta_n;
	double MS = ef.G[ind].M0 + nn * t;
	double E = MS;
	double rez = E + 1;
	while (fabs(rez - E) > 1e-10) {
		rez = E;
		E = MS + ef.G[ind].e * sin(E);
	}
	double v = atan2(SQRT(1 - SQR(ef.G[ind].e)) * sin(E), cos(E) - ef.G[ind].e);
	double f = v + ef.G[ind].w0;
	double delu = ef.G[ind].Cus * sin(2 * f) + ef.G[ind].Cuc * cos(2 * f);
	double delr = ef.G[ind].Crs * sin(2 * f) + ef.G[ind].Crc * cos(2 * f);
	double deli = ef.G[ind].Cis * sin(2 * f) + ef.G[ind].Cic * cos(2 * f);
	double u = f + delu;
	double r = SQR(ef.G[ind].A) * (1 - ef.G[ind].e * cos(E)) + delr;
	double ii = ef.G[ind].i0 + ef.G[ind].IDOT * t + deli;
	double ugl = ef.G[ind].Omega0 + t * (ef.G[ind].Omega_p - UGLE) - UGLE * ef.G[ind].Toe;
	
	double Xorb = r * cos(u);
	double Yorb = r * sin(u);

	double x = Xorb * cos(ugl) - Yorb * cos(ii) * sin(ugl);
	double y = Xorb * sin(ugl) + Yorb * cos(ii) * cos(ugl);
	double z = Yorb * sin(ii);
	double dT = 0;

	double Toc = time2gpst(epoch2time(ef.time), NULL);
	dT += -2. * SQRT(MU) * ef.G[ind].A * ef.G[ind].e * sin(E) / SQR(SPEED_OF_LIGHT);//релятивистский эффект
	dT += ef.G[ind].af0;//отклонение часов спутника
	dT += ef.G[ind].af1 * (Tem - Toc);
	dT += ef.G[ind].af2 * (Tem - Toc) * (Tem - Toc);

	return (SatParams) { .X = x, .Y = y, .Z = z, .dT = dT };
}

//частота в зависимости от типа спутника и канала частоты 
inline double f1(char* typeOfSat, int K) { return(typeOfSat == "G") ? GF1 : RF1 + K * DRF1; }

inline double f2(char* typeOfSat, int K) { return(typeOfSat == "G") ? GF2 : RF2 + K * DRF2; }

SatParams RSatPosition(GFile ef, int ind, Date date) {
	double q[6] = { ef.R[ind].X,ef.R[ind].Y,ef.R[ind].Z,ef.R[ind].Vx,ef.R[ind].Vy, ef.R[ind].Vz };

	double efsec = epoch2time(ef.time) + ef.R[ind].leap;
	double dateSec= epoch2time(date);

	RK4(model, q, q, efsec, dateSec, 10, 6);

	double dT = 0;
	dT += ef.R[ind].TauN;
	dT += ef.R[ind].TauC;
	dT += ef.R[ind].GammaN * (dateSec - efsec);

	return (SatParams){ .X = q[0],.Y = q[1], .Z = q[2], .dT = dT };
}

SatParams SatPosition(navFile ef, const char* typeOfSat, int ind, Date date) {
	errexit(typeOfSat != "G" && typeOfSat != "R", "ERROR(SatPosition):incorrect type of satellite");
	return typeOfSat == "G" ? GSatPosition(ef.n, ind, date) : RSatPosition(ef.g, ind, date);
}

void satPositions(navFiles* F, navFiles* prev,char*typeOfsat, SatVector* S) {
	errexit(typeOfsat != "G" && typeOfsat != "R", "ERROR(satPositions):incorrect type of satellite");

	int NUM = (typeOfsat == "G") ? GNUM : RNUM;

	for (int i = 0; i < F->o.size; i++) {
		Sat satList;
		int j;
		Date _date = F->o.container[i].time;

		
		//распараллеливание независимых итераций
		#pragma omp parallel for private(j)
		for (j = 1; j <= NUM; j++) {
			OFileParams SatFromO;
			satList.sat[j].satVisible = false;
			int prevSize, navSize;
			if (typeOfsat == "R") {
				prevSize = prev->g.size;
				navSize = F->g.size;
				SatFromO = F->o.container[i].sat.R[j];
			}
			else {
				prevSize = prev->n.size;
				navSize = F->n.size;
				SatFromO = F->o.container[i].sat.G[j];
			}
			if (SatFromO.satVisible) {
				for (int k = prevSize + navSize - 1; k > -1; k--) {
					navFile _lastEf;
					bool navSatVisible;
					Date lastEfTime;
					if (typeOfsat == "R") {
						_lastEf.g= (k < prevSize) ? prev->g.container[k] : F->g.container[k - prevSize];
						navSatVisible = _lastEf.g.R[j].satVisible;
						lastEfTime = _lastEf.g.time;
					}
					else {
						_lastEf.n = (k < prevSize) ? prev->n.container[k] : F->n.container[k - prevSize];
						navSatVisible = _lastEf.n.G[j].satVisible;
						lastEfTime = _lastEf.n.time;
					}

					double lastefsec, oSec;
					lastefsec = epoch2time(lastEfTime);
					oSec = epoch2time(_date);

					if (navSatVisible && (lastefsec <= oSec)) {
						Date trDate = _date;

						trDate.seconds -= SatFromO.p1 / SPEED_OF_LIGHT;
						//рассчёт координат навигационного спутника
						satList.sat[j] = SatPosition(_lastEf, typeOfsat, j, trDate);

						//добавление дополнительных параметров для дальнейшего определение координат станции и вычисления ПЭС
						satList.sat[j].satVisible = true;
						satList.sat[j].P1 = SatFromO.p1;
						satList.sat[j].P2 = SatFromO.p2;
						satList.sat[j].L1 = SatFromO.l1;
						satList.sat[j].L2 = SatFromO.l2;
						satList.sat[j].C1 = SatFromO.c1;
						satList.sat[j].C2 = SatFromO.c2;

						satList.sat[j].K = (typeOfsat == "R") ? _lastEf.g.R[j].K : 0;
						double f1_ = f1(typeOfsat, satList.sat[j].K);
						double f2_ = f2(typeOfsat, satList.sat[j].K);
						satList.sat[j].P = (SatFromO.p1 * SQR(f1_) - SatFromO.p2 * SQR(f2_)) / (SQR(f1_) - SQR(f2_));
						satList.sat[j].C = (SatFromO.c1 * SQR(f1_) - SatFromO.c2 * SQR(f2_)) / (SQR(f1_) - SQR(f2_));
						
						//satList.sat[j].P = SatFromO.p1;

						break;
					}
				}
			}
		}
		satList.time = _date;
		SatPush(S, &satList);
	}
	if (typeOfsat == "R") GFileCopy(&(prev->g), &(F->g));
	else if (typeOfsat == "G") NFileCopy(&(prev->n), &(F->n));
}



int pointPosition(Sat sat, char* typeOfSat, double* pointPos) {
	errexit(typeOfSat != "G" && typeOfSat != "R", "ERROR(pointPosition):incorrect type of satellite");

	int N = (typeOfSat == "G") ? GNUM : RNUM;

	SatParamsVector visibleSat;
	SatParamsInit(&visibleSat, 1);
	//выделение тех спутников,которые видны в текущий момент времени 
	for (int j = 1; j <= N; j++) {
		SatParams _sat;
		if (sat.sat[j].satVisible) {
			_sat = sat.sat[j];
			SatParamsPush(&visibleSat, &_sat);
		}
	}
		//если видимых спутников меньше 4, то определить местоположение невозможно
	if (visibleSat.size < 4) return 1;


	Matrix* X0 = MatrixCreate(4, 1);

	MatrixSet(X0, 0, 0, 0);
	MatrixSet(X0, 1, 0, 0);
	MatrixSet(X0, 2, 0, 0);
	MatrixSet(X0, 3, 0, 0);
	int lsqiter = 0;

	while (lsqiter < LSQITERMAX) {
		Matrix* A = MatrixCreate(visibleSat.size, 4);
		Matrix* B = MatrixCreate(visibleSat.size, 1);
		lsqiter++;

		for (int row = 0; row < visibleSat.size; row++) {

			double RsRr[3] = {
				visibleSat.container[row].X - MatrixGetV(X0, 0, 0),
				visibleSat.container[row].Y - MatrixGetV(X0, 1, 0),
				visibleSat.container[row].Z - MatrixGetV(X0, 2, 0)
			};

			double roRsRr = l2norm(RsRr, 3);

			for(int i=0;i<3;i++)
				MatrixSet(A, row, i, -RsRr[i] / roRsRr);
			MatrixSet(A, row, 3, 1);

			roRsRr += UGLE * (visibleSat.container[row].X * MatrixGetV(X0, 1, 0) - visibleSat.container[row].Y * MatrixGetV(X0, 0, 0)) / SPEED_OF_LIGHT;

			roRsRr -= visibleSat.container[row].dT * SPEED_OF_LIGHT;

			roRsRr += MatrixGetV(X0, 3, 0);

			MatrixSet(B, row, 0, visibleSat.container[row].P - roRsRr);
		}

		Matrix* deltaX = lsq(A, B);

		Matrix* X1 = MatrixSum(X0, deltaX);
		MatrixFree(X0); free(X0);
		X0 = X1;

		double delta[4] = {
			MatrixGetV(deltaX, 0, 0),
			MatrixGetV(deltaX, 1, 0),
			MatrixGetV(deltaX, 2, 0),
			MatrixGetV(deltaX, 3, 0)
		};

		MatrixFree(A); free(A);
		MatrixFree(B); free(B);
		MatrixFree(deltaX); free(deltaX);

		if (linfnorm(delta, 4) < 1e-5) break;
	}

	errexit(lsqiter >= LSQITERMAX, "ERROR(pointPosition):lsq iteration owerflow");

	pointPos[0] = MatrixGetV(X0, 0, 0);
	pointPos[1] = MatrixGetV(X0, 1, 0);
	pointPos[2] = MatrixGetV(X0, 2, 0);

	if (fabs(MatrixGetV(X0, 3, 0)/SPEED_OF_LIGHT) > 0.005) return 1;

	MatrixFree(X0);
	free(X0);
	SatParamsFree(&visibleSat);
	
	return 0;
}

void TEC(SatVector *sat,double *pointPos, char * typeOfSat,FILE *All,FILE **nav) {
	errexit(typeOfSat != "G" && typeOfSat != "R", "ERROR(pointPosition):incorrect type of satellite");

	int Num= typeOfSat == "G"? GNUM: RNUM;

	double Rr = l2norm(pointPos, 3);

	double tecPN = 0, tecLN =0;
	double NN = 0;

	for (int i = 0; i < sat->size; i++) {
		double tecC = 0;
		double tecL = 0;

		int N = 0;

		for (int j = 1; j <= Num; j++) {
			if (sat->container[i].sat[j].satVisible) {
				double RsRr[3] = {
					sat->container[i].sat[j].X - pointPos[0],
					sat->container[i].sat[j].Y - pointPos[1],
					sat->container[i].sat[j].Z - pointPos[2]
				};
				double cosz = (RsRr[0]*pointPos[0] + RsRr[1] *pointPos[1] + RsRr[2] *pointPos[2]) / (Rr*l2norm(RsRr,3));
				if (asin(cosz) * 180 / PI < 10) continue;
				int K = sat->container[i].sat[j].K;
				double cosz_ = pow((1 - (1 - cosz * cosz) * pow(Rr / (RC + H_ION), 2)), 0.5);

				double beta = SQR(f1(typeOfSat, K) * f2(typeOfSat, K)) / (K_ION * 1e+16 * (SQR(f1(typeOfSat, K)) - SQR(f2(typeOfSat, K))));

				double _tecC = beta * (sat->container[i].sat[j].P2 - sat->container[i].sat[j].P1) * cosz_;
				
				double _tecL = beta * (sat->container[i].sat[j].L1 * SPEED_OF_LIGHT / f1(typeOfSat, K) - sat->container[i].sat[j].L2 * SPEED_OF_LIGHT / f2(typeOfSat, K)) * cosz_;

				if (typeOfSat == "R") _tecL *= -1;


				fprintf(nav[j - 1], "%d/%d/%d %d:%d:%d//%15lf//%15lf//%f//%f\n", sat->container[i].time.day, sat->container[i].time.month, sat->container[i].time.year, sat->container[i].time.hours, sat->container[i].time.minutes, sat->container[i].time.seconds, _tecC, _tecL, sat->container[i].sat[j].P1 - sat->container[i].sat[j].P, sat->container[i].sat[j].P2 - sat->container[i].sat[j].P);

				tecC += _tecC;
				tecL += _tecL;

				N++;
			}
		}
		
		if (N > 0) {
			tecC /= N;
			tecL /= N;
			tecPN += tecC;
			tecLN += tecL;
			NN++;
			if (NN > 0) {
				tecPN /= NN;
				tecLN /= NN;
				fprintf(All, "%d/%d/%d %d:%d:%d//%15lf//%15lf\n", sat->container[i].time.day, sat->container[i].time.month, sat->container[i].time.year, sat->container[i].time.hours, sat->container[i].time.minutes, sat->container[i].time.seconds, tecPN, tecLN);
				NN = 0;
				tecPN = 0;
				tecLN = 0;
			}
		
		}
	}
}

//вывод карты ПЭС в файл, аналогично поиску ПЭС по измерениям ищем ПЭС для всех сптников в области видимости на каждый момент времени
void modelTEC(SatVector *sat, double *pointPos, char * typeOfSat, IonexFileVector *iF, FILE *Map){
	errexit(typeOfSat != "G" && typeOfSat != "R", "ERROR(modelTEC):incorrect type of satellite");

	int Num = typeOfSat == "G" ? GNUM : RNUM;

	for (int i = 0; i < sat->size; i++) {
		double tec = 0;

		int N = 0;
		for (int j = 1; j <= Num; j++) {
			if (sat->container[i].sat[j].satVisible) {
				double func = ((sat->container[i].sat[j].X - pointPos[0])*pointPos[0] + (sat->container[i].sat[j].Y - pointPos[1])*pointPos[1] + (sat->container[i].sat[j].Z - pointPos[2])*pointPos[2]) / (fabs(sqrt(pointPos[0] * pointPos[0] + pointPos[1] * pointPos[1] + pointPos[2] * pointPos[2])*sqrt(pow((sat->container[i].sat[j].X - pointPos[0]), 2) + pow((sat->container[i].sat[j].Y - pointPos[1]), 2) + pow((sat->container[i].sat[j].Z - pointPos[2]), 2))));
				if (asin(func) * 180 / PI < 10) continue;//если КА ниже 10 градусов не учитываем его
				//находим точку пересечения прямой соединяющую станцию и КА со сферой ПЭС, для этого решаем квадратное уравнение A*t^2+B*t+C=0
				double A, B, C;
				A = pow(sat->container[i].sat[j].X - pointPos[0], 2) + pow(sat->container[i].sat[j].Y - pointPos[1], 2) + pow(sat->container[i].sat[j].Z - pointPos[2], 2);
				B = 2 * ((sat->container[i].sat[j].X - pointPos[0])*pointPos[0] + (sat->container[i].sat[j].Y - pointPos[1])*pointPos[1] + (sat->container[i].sat[j].Z - pointPos[2])*pointPos[2]);
				C = pow(pointPos[0], 2) + pow(pointPos[1], 2) + pow(pointPos[2], 2) - pow(RC + H_ION, 2);
				double t[2];
				t[0] = (-B + sqrt(B*B - 4 * A*C)) / (2 * A);
				t[1] = (-B - sqrt(B*B - 4 * A*C)) / (2 * A);
				double Xtec[3], Ytec[3], Ztec[3];
				//нахоим две точки пересечения 
				Xtec[0] = (sat->container[i].sat[j].X - pointPos[0])*t[0] + pointPos[0];
				Ytec[0] = (sat->container[i].sat[j].Y - pointPos[1])*t[0] + pointPos[1];
				Ztec[0] = (sat->container[i].sat[j].Z - pointPos[2])*t[0] + pointPos[2];
				Xtec[1] = (sat->container[i].sat[j].X - pointPos[0])*t[1] + pointPos[0];
				Ytec[1] = (sat->container[i].sat[j].Y - pointPos[1])*t[1] + pointPos[1];
				Ztec[1] = (sat->container[i].sat[j].Z - pointPos[2])*t[1] + pointPos[2];
				//отбрасываем ту, которая находится дальше, потому что из геометрического смысла нам подходит та, что ближе
				if (sqrt(pow(Xtec[0] - pointPos[0], 2) + pow(Ytec[0] - pointPos[1], 2) + pow(Ztec[0] - pointPos[2], 2)) < sqrt(pow(Xtec[1] - pointPos[0], 2) + pow(Ytec[1] - pointPos[1], 2) + pow(Ztec[1] - pointPos[2], 2))) {
					Xtec[2] = Xtec[0];
					Ytec[2] = Ytec[0];
					Ztec[2] = Ztec[0];
				}
				else {
					Xtec[2] = Xtec[1];
					Ytec[2] = Ytec[1];
					Ztec[2] = Ztec[1];
				}
				double tecm = TECmapDec(iF, sat->container[i].time, Xtec[2], Ytec[2], Ztec[2]);//рассчитываем ПЭС в точке пересечения ионосферы в момент когда происходило измерение псевдодальности до выбранного КА
				if (tecm != 9999) {// если 9999 значит точка выходит за пределы карты, либо нет данных
					tec += tecm;
					N++;
				}
			}

		}
		if (N > 0) {
			tec /= N;

			fprintf(Map, "%d/%d/%d %d:%d:%d//%15lf\n", sat->container[i].time.day, sat->container[i].time.month, sat->container[i].time.year, sat->container[i].time.hours, sat->container[i].time.minutes, sat->container[i].time.seconds, tec);
		}
	}
}


int main() {

	//определение числа потоков при параллельных вычислениях 
	omp_set_num_threads(3);

	StrVector oFilenameList, nFilenameList, gFilenameList, codgFilenameList;
	
	StrInit(&oFilenameList, 365, 20);
	StrInit(&nFilenameList, 365, 20);
	StrInit(&gFilenameList, 365, 20);
	StrInit(&codgFilenameList, 365, 20);

	const char *oPattern = "*.*o";
	const char *nPattern = "*.*n";
	const char *gPattern = "*.*g";
	const char *mapPattern = "codg*.*i";
	
	findAllFilename(&oFilenameList, oPattern);
	findAllFilename(&nFilenameList, nPattern);
	findAllFilename(&gFilenameList, gPattern);
	findAllFilename(&codgFilenameList, mapPattern);
	
	const char Gfiles[GNUM][8] = { "G01.txt","G02.txt","G03.txt","G04.txt","G05.txt","G06.txt","G07.txt","G08.txt","G09.txt","G10.txt",
		"G11.txt","G12.txt","G13.txt","G14.txt","G15.txt","G16.txt","G17.txt","G18.txt","G19.txt","G20.txt",
		"G21.txt","G22.txt","G23.txt","G24.txt","G25.txt","G26.txt","G27.txt","G28.txt","G29.txt","G30.txt",
		"G31.txt","G32.txt" };
	const char Rfiles[RNUM][8] = { "R01.txt","R02.txt","R03.txt","R04.txt","R05.txt","R06.txt","R07.txt","R08.txt","R09.txt","R10.txt",
		"R11.txt","R12.txt","R13.txt","R14.txt","R15.txt","R16.txt","R17.txt","R18.txt","R19.txt","R20.txt",
		"R21.txt","R22.txt","R23.txt","R24.txt" };

	FILE *allG;
	allG = fopen("allG.txt", "w");

	FILE *allR;
	allR = fopen("allR.txt", "w");

	FILE *TECMapG;
	TECMapG = fopen("TECMapG.txt", "w");

	FILE *TECMapR;
	TECMapR = fopen("TECMapR.txt", "w");
	FILE* pos;
	pos = fopen("pos.txt", "w");

	FILE *GTEC[GNUM];
	for (int i = 0; i < GNUM; i++) GTEC[i] = fopen(Gfiles[i], "w");
	
	FILE *RTEC[RNUM];
	for (int i = 0; i < RNUM; i++) RTEC[i] = fopen(Rfiles[i], "w");
	
	navFiles prev;

	NFileInit(&prev.n, 0);
	GFileInit(&prev.g, 0);

	for (int ofileNumber = 0; ofileNumber < 5; ofileNumber++) {

		Date oDate = getFileDate(oFilenameList.container[ofileNumber]);

		int nfileNumber=0, gfileNumber=0;

		for (; nfileNumber < nFilenameList.size; nfileNumber++) {
			Date nDate = getFileDate(nFilenameList.container[nfileNumber]);
			if (cmpDate(&nDate, &oDate)) break;
		}
		for (; gfileNumber < gFilenameList.size; gfileNumber++) {
			Date gDate = getFileDate(gFilenameList.container[gfileNumber]);
			if (cmpDate(&gDate, &oDate)) break;
		}

		if (nfileNumber == nFilenameList.size) {
			NFileFree(&prev.n);
			NFileInit(&prev.n, 0);
			printf("%i/%i/%i GPS ephemeris data not found\n", oDate.year, oDate.month, oDate.day);
			continue;
		}
		if (gfileNumber == gFilenameList.size) {
			GFileFree(&prev.g);
			GFileInit(&prev.g, 0);
			printf("%i/%i/%i GLONASS ephemeris data not found\n", oDate.year, oDate.month, oDate.day);
			continue;
		}

		long t = clock();
		navFiles F;
		IonexFileVector iF;
		SatVector G,R;
		double pointPos[3];

		SatInit(&G, 2881);
		SatInit(&R, 2881);

		//распараллеливание независимых блоков программы, пасинг файлов
		#pragma omp parallel sections num_threads(4)
		{
			#pragma omp section
			{
				OFileInit(&F.o, 2881);
				parseOFile(&oFilenameList, ofileNumber, &F.o);
			}
			#pragma omp section
			{		
				NFileInit(&F.n, 14);
				parseNFile(&nFilenameList, nfileNumber, &F.n);
			}
			#pragma omp section
			{
				GFileInit(&F.g, 48);
				parseGFile(&gFilenameList, gfileNumber, &F.g);
			}
			#pragma omp section
			{
				IonexFileInit(&iF, 25);
				parseIonex(&codgFilenameList, ofileNumber, &iF);
			}
		}
		double _year = F.o.container[0].time.year;
		double _month = F.o.container[0].time.month;
		double _day = F.o.container[0].time.day;

		errexit(iF.container[0].time.year != _year || iF.container[0].time.month != _month || iF.container[0].time.day != _day, "IONEX data not found");
		
		satPositions(&F, &prev, "G", &G);
		satPositions(&F, &prev, "R", &R);

		/**/
		for (int i = 0; i < G.size + R.size; i++)
			if (i < G.size) {
				if (pointPosition(G.container[i], "G", pointPos) == 0) break;
			}
			else {
				if (pointPosition(R.container[i - G.size], "R", pointPos) == 0) break;
			}
		/**/
		//распараллеливание независимых блоков программы, определение ПЭС для ГЛОНАСС и GPS
		#pragma omp parallel sections num_threads(4)
		{
			#pragma omp section
			{
				TEC(&G, pointPos, "G", allG, GTEC);
			}
			#pragma omp section
			{
				TEC(&R, pointPos, "R", allR, RTEC);
			}
			#pragma omp section
			{
				modelTEC(&G, pointPos, "G", &iF, TECMapG);
			}
			#pragma omp section
			{
				modelTEC(&R, pointPos, "R", &iF, TECMapR);
			}
		}
	
		SatFree(&G);SatFree(&R);
		OFileFree(&F.o); NFileFree(&F.n); GFileFree(&F.g); 
		IonexFileFree(&iF);
		t = clock() - t;
		printf("%i\n", t);
	}

	StrFree(&oFilenameList); StrFree(&nFilenameList); StrFree(&gFilenameList);
	NFileFree(&prev.n); GFileFree(&prev.g);

	fclose(allG); fclose(allR); fclose(TECMapG); fclose(TECMapR);
	for (int i = 0; i < GNUM; i++) fclose(GTEC[i]);
	for (int i = 0; i < RNUM; i++) fclose(RTEC[i]);

	return 0;
}