#include "LatLonToUTM.h"
#include <math.h>

namespace LL2UTM{
	const Ellipse ellipse[3] = { 
		Ellipse(6378137.0, 298.257223563, "WGS84"), 
		Ellipse(6378137.0, 298.257222101, "GRS80"),
		Ellipse(6377397.155, 299.1528128, "Bessel1841")
	};

	void GeoConverter::selectEllipseModel(EllipseModel model){
		ellipse_model = model;
		double f = 1.0 / ellipse[ellipse_model].f_inv;
		e2 = 2 * f - f*f;
		e4 = e2 * e2;
		e6 = e4 * e2;
		he2 = e2 / (1.0 - e2);
	}

	const Ellipse& GeoConverter::getEllipseModel() const{
		return ellipse[ellipse_model];
	}

	bool GeoConverter::findZone(double lat, double lon, double& refLon, UTMZone& zone){
		const static char LAT_ZONE_CHAR[] = "CDEFGHJKLMNPQRSTUVWX";//20 zones along with lat
		double lon_tmp = lon; 
		double lat_tmp = lat;
		bool exceptional_zone = false;

#pragma region Invalid_Lat_Lon
		if (fabs(lon_tmp) > 180.0 || fabs(lat_tmp) > 90.0){
			refLon = -1.0;
			zone = UTMZone();
			return false;
		}
#pragma endregion


#pragma region Antarctic_And_Arctic_Zone
		if (lat_tmp >= 84.0){
			if (lon_tmp >= 0.0){
				zone.latDiv = 'Z';
				zone.lonDiv = -1;
				refLon = 0.0;
			}
			else{
				zone.latDiv = 'Y';
				zone.lonDiv = -1;
				refLon = 180.0;
			}
			return true;
		}
		else if (lat_tmp < -80.0){
			if (lon_tmp >= 0.0){
				zone.latDiv = 'B';
				zone.lonDiv = -1;
				refLon = 0.0;
			}
			else{
				zone.latDiv = 'A';
				zone.lonDiv = -1;
				refLon = 180.0;
			}
			return true;
		}
#pragma endregion


		if (lat_tmp >= 80.0) zone.latDiv = LAT_ZONE_CHAR[19];//Zone X ranges [72, 84)
		else                 zone.latDiv = LAT_ZONE_CHAR[int((lat_tmp + 80) / 8)];

		zone.lonDiv = int((lon_tmp + 180) / 6) + 1; 
		
		if (zone.lonDiv == 61){
			//When lon == 180.0(exactly) and EW == 'E'
			zone.lonDiv = 60;
		}
#pragma region Exceptional_Zone_Cases
		if (zone.latDiv == 'V'){
			if (lon_tmp >= 0.0 && lon_tmp < 12.0){
				exceptional_zone = true;
				if (lon_tmp < 3.0){
					refLon = 1.5;
				}
				else{
					//Exception at Norway's southwest coast(zone 32V) -> [3.0~12.0) degrees
					zone.lonDiv = 32;
					refLon = 7.5;
				}
			}
		}
		if (zone.latDiv == 'X'){
			//Exceptions at Svalbard
			if (lon_tmp >= 0.0 && lon_tmp < 42.0){
				exceptional_zone = true;
				if (lon_tmp < 9.0) {
					zone.lonDiv = 31; //31X
					refLon = 4.5;
				}
				else if (lon_tmp < 21.0) {
					zone.lonDiv = 33; //32X
					refLon = 15.0;
				}
				else if (lon_tmp < 33.0) {
					zone.lonDiv = 35; //35X
					refLon = 27.0;
				}
				else {
					zone.lonDiv = 37; //37X
					refLon = 37.5;
				}
			}
		}
#pragma endregion

		if (!exceptional_zone){
			if (zone.lonDiv > 30){
				refLon = (zone.lonDiv - 31) * 6 + 3;
			}
			else{
				refLon = (zone.lonDiv - 31) * 6 + 3;
			}
		}
		return true;
	}

	void GeoConverter::convertLL2UTM(double lat, double lon, UTMZone& zone, double& easting, double& northing, double& k){
#pragma region Saskatchewan_Univ_Version
		//Equation from (http://www.engr.usask.ca/classes/CE/316/notes/CE%20316%20CH%204C%2031-1-12%20-INSTRUCTOR.pdf)
		double k0 = 0.9996;
		double ref_lon;
		findZone(lat, lon, ref_lon, zone);
		lat = deg2rad(lat);
		lon = deg2rad(lon);
		ref_lon = deg2rad(ref_lon);
		
		// Start calcuation.
		const double N0 = lat < 0.0 ? 10000000 : 0.0;
		const double E0 = 500000.0;
		double sin_lat = sin(lat);
		double sin_lon = sin(lon);
		double cos_lon = cos(lon);
		double cos_lat = cos(lat);
		double Rm = ellipse[ellipse_model].a * (1.0 - e2) / pow(1 - e2*sin_lat*sin_lat, 3.0 / 2.0);
		double Rn = ellipse[ellipse_model].a / sqrt(1 - e2*sin_lat*sin_lat);
		double T = tan(lat)*tan(lat);
		double C = he2*cos_lat*cos_lat;
		double A = (lon - ref_lon) * cos_lat;
		double AA = A*A;
		double AAA = AA*A;
		double M = ellipse[ellipse_model].a * ((1 - e2 / 4 - 3 * e4 / 64 - 5 * e6 / 256)*lat - (3.0 / 8 * e2 + 3.0 / 32 * e4 + 45.0 / 1024 * e6)*sin(2 * lat)
			+ (15.0 / 256 * e4 + 45.0 / 1024 * e6)*sin(4 * lat) - (35.0 / 3072 * e6)*sin(6 * lat));
		easting = E0 + k0 * Rn *(A + (1 - T + C)*AAA / 6 + (5 - 18.0*T + T*T + 72 * C - 58 * he2)*AAA*AA / 120.0);
		northing = N0 + k0 * (M + Rn*tan(lat)*(AA / 2 + (5 - T + 9 * C + 4 * C*C)*AA*AA / 24 + (61.0 - 58 * T + T*T + 600 * C - 330 * he2)*AAA*AAA / 720.0));
		k = k0 * (1.0 + (1 + C)*AA / 2 + (5 - 4 * T + 42 * C + 13 * C*C - 28 * he2)*AA*AA / 24 + (61 - 148 * T + 16 * T*T)*AAA*AAA / 720);
#pragma endregion


#pragma region Wikipedia_Version
		/*
//Equation from Wikipedia(https://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system)
//Not working.. Why?
double n = 1 / (ellipse[ellipse_model].f_inv*2.0 - 1.0);
double nn = n*n;
double nnn = nn*n;
double A = ellipse[ellipse_model].a / (1 + n) * (1.0 + nn / 4.0 + nn*nn / 64.0 + nn*nn*nn / 256.0);
double a1 = 1.0 / 2 * n - 2.0 / 3 * nn + 5.0 / 16 * nnn;
double a2 = 13.0 / 48 * nn - 3.0 / 5 * nnn;
double a3 = 61.0 / 240 * nnn;
double b1 = 1.0 / 2 * n - 2.0 / 3 * nn + 37.0 / 96 * nnn;
double b2 = 1.0 / 48 * nn + 1.0 / 15 * nnn;
double b3 = 17.0 / 480 * nnn;
double g1 = 2.0 * n - 2.0 / 3 * nn - 2 * nnn;
double g2 = 7.0 / 3 * nn - 8.0 / 5 * nnn;
double g3 = 56.0 / 15 * nnn;
double ref_lon;

findZone(lat, lon, ref_lon, zone);

lat = deg2rad(lat);
lon = deg2rad(lon);
ref_lon = deg2rad(ref_lon);
double t = sinh(atanh(sin(lat)) - 2 * sqrt(n) / (1 + n)*atanh(2 * sqrt(n) / (1 + n)*sin(lat)));
double ksai = atan(t / cos(t/(lon - ref_lon)));
double eta = atanh(sin(lon - ref_lon) / sqrt(1 + t*t));
double sigma = 1 + 2 * (a1*sin(2 * ksai)*cosh(2 * eta) + 2 * a2*sin(4 * ksai)*cosh(4 * eta) + 3 * a3*sin(6 * ksai)*cosh(6 * eta));
double tau = 2*(a1*sin(2*ksai)*sinh(2*eta) + 2*a2*sin(4*ksai)*sinh(4*eta) + 3*a3*sin(6*ksai)*sinh(6*eta));

easting = 500 + 0.9996*A*(eta + a1*cos(2 * ksai)*sinh(2 * eta) + a2*cos(4 * ksai)*sinh(4 * eta) + a3*cos(6 * ksai)*sinh(6 * eta));
northing = 0 + 0.9996*A*(ksai + a1*sin(2 * ksai)*cosh(2 * eta) + a2*sin(4 * ksai)*cosh(4 * eta) + a3*sin(6 * ksai)*cosh(6 * eta));
*/
#pragma endregion

	}
}