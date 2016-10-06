// Copyright (c) 2016 JaeJunLee
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
// DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
// DEALINGS IN THE SOFTWARE.
#include <gtest\gtest.h>
#define private public
#define protected public
#include "..\UTMCoordinates\LatLonToUTM.h"
#include <math.h>

#define ZONE_TEST_CASE(in_lat, in_lon, ref_lon_des, lat_zone_des, lon_zone_des) do{\
		LL2UTM::GeoConverter::findZone(in_lat, in_lon, ref_lon, zone); \
		EXPECT_DOUBLE_EQ(ref_lon_des, ref_lon); \
		EXPECT_EQ(lat_zone_des, zone.latDiv); \
		EXPECT_EQ(lon_zone_des, zone.lonDiv); \
	}while (0)

TEST(GeoConverter, findZone){
	LL2UTM::UTMZone zone;
	double ref_lon;
	
	//ZONE_TEST_CASE(lat, lon, ref_lon_des, lat_zone_des, lon_zone_des);
	//North - East
	ZONE_TEST_CASE(1.0, 1.0,  3.0, 'N', 31);
	ZONE_TEST_CASE(0.001, 0.001,  3.0, 'N', 31);
	ZONE_TEST_CASE(35.9078, 127.7669,  129, 'S', 52);
	ZONE_TEST_CASE(56.0, 1.0, 1.5, 'V', 31);//31V
	ZONE_TEST_CASE(63.0, 6.0,  7.5, 'V', 32);//32V
	ZONE_TEST_CASE(60.0, 11.9,  7.5, 'V', 32);//32V
	ZONE_TEST_CASE(72.0, 1.0,  4.5, 'X', 31);//31X
	ZONE_TEST_CASE(75.0, 9.5,  15.0, 'X', 33);//33X
	ZONE_TEST_CASE(80.0, 24.0,  27.0, 'X', 35);//35X
	ZONE_TEST_CASE(83.9, 34.0,  37.5, 'X', 37);//37X
	
	//
	ZONE_TEST_CASE(-0.001, -0.001, -3.0, 'M', 30);
	ZONE_TEST_CASE(-1.0,  -1.0, -3.0, 'M', 30);
	ZONE_TEST_CASE(-15.0, -7.0, -9.0, 'L', 29);
}

/*
	Rounding-ups are performed to expected and actual at 0.1^(n_prec+1).
	And if the difference between two rounded values is zero, they are considered same.

	Be warn that doubleEQWithTolerance(1.1532, 1.1536, 3) will return false, 
	1.1532 --> 1.1530
	1.1536 --> 1.1540
*/
bool doubleEQWithTolerance(double expected, double actual, const int& n_prec){
	double mult = pow(10.0, n_prec);
	int int_expected = int(expected * mult);
	int int_actual = int(actual * mult);
	if (int_expected-int_actual == 0){
		return true;
	}
	return false;
}

TEST(Auxilary, doubleEQWithTolerance){
	EXPECT_TRUE(doubleEQWithTolerance(1.1532, 1.1531, 1));
	EXPECT_TRUE(doubleEQWithTolerance(1.1532, 1.1531, 2));
	EXPECT_TRUE(doubleEQWithTolerance(1.1532, 1.1531, 3));
	EXPECT_FALSE(doubleEQWithTolerance(1.1532, 1.1531, 4));

	EXPECT_TRUE(doubleEQWithTolerance(1.1532, 1.1536, 1));
	EXPECT_TRUE(doubleEQWithTolerance(1.1532, 1.1536, 2));
	EXPECT_TRUE(doubleEQWithTolerance(1.1532, 1.1536, 3));//Be warn that this will return false, because of round-up
	EXPECT_FALSE(doubleEQWithTolerance(1.1532, 1.1536, 4));
}
#if 1
	//flooring
	#define LL2UTM_TEST_CASE(in_lat, in_lon, easting_des, northing_des, n_prec) do{\
		converter.convertLL2UTM(in_lat, in_lon, zone, easting, northing, k); \
		double mult = pow(10, n_prec);\
		EXPECT_DOUBLE_EQ(floor(easting_des*mult)/mult, floor(easting*mult)/mult); \
		EXPECT_DOUBLE_EQ(floor(northing_des*mult)/mult, floor(northing*mult)/mult); \
	}while (0)
#else
	//rounding
	#define LL2UTM_TEST_CASE(in_lat, in_lon, easting_des, northing_des, n_prec) do{\
		converter.convertLL2UTM(in_lat, in_lon, zone, easting, northing, k); \
		double mult = pow(10, n_prec); \
		EXPECT_DOUBLE_EQ(round(easting_des*mult) / mult, round(easting*mult) / mult); \
		EXPECT_DOUBLE_EQ(round(northing_des*mult) / mult, round(northing*mult) / mult); \
	}while (0)
#endif
TEST(GeoConverter_convertLL2UTM, NorthEast){
	LL2UTM::GeoConverter converter;
	LL2UTM::UTMZone zone;
	converter.selectEllipseModel(LL2UTM::GeoConverter::WGS84);
	double easting, northing, k;
	LL2UTM_TEST_CASE(35.9078, 127.7669, 388731.7, 3974424.5, 0); //TODO : we should find accurate reference for Korea.
	LL2UTM_TEST_CASE(1.0, 1.0, 277438.2635, 110597.9725, 3);
	LL2UTM_TEST_CASE(3.1548, 15.1548, 517199.4544, 348705.0046,3);
}

TEST(GeoConverter_convertLL2UTM, NorthWest){
	LL2UTM::GeoConverter converter;
	LL2UTM::UTMZone zone;
	converter.selectEllipseModel(LL2UTM::GeoConverter::WGS84);
	double easting, northing, k;
	LL2UTM_TEST_CASE(10.0, -3.5, 445201.5878, 1105454.0118, 3);
}

TEST(GeoConverter_convertLL2UTM, SouthWest){
	LL2UTM::GeoConverter converter;
	LL2UTM::UTMZone zone;
	converter.selectEllipseModel(LL2UTM::GeoConverter::WGS84);
	double easting, northing, k;
	LL2UTM_TEST_CASE(-15.0, -3.5, 446245.6708, 8341613.3000, 3);
}

TEST(GeoConverter_convertLL2UTM, SouthEast){
	LL2UTM::GeoConverter converter;
	LL2UTM::UTMZone zone;
	converter.selectEllipseModel(LL2UTM::GeoConverter::WGS84);
	double easting, northing, k;
	LL2UTM_TEST_CASE(-15.0, 4.5, 661277.2660, 8341127.5191, 3);
	
}

int main(int argc, char** argv){
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}