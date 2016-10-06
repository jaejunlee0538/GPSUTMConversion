// Copyright (c) 2016 JaeJunLee
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
// DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
// DEALINGS IN THE SOFTWARE.

#ifndef LatLonToUTM_H_
#define LatLonToUTM_H_
#include <stdint.h>
#include <string>

namespace LL2UTM{
	struct UTMZone{
		char latDiv;
		int8_t lonDiv;
		UTMZone():latDiv(0),lonDiv(-1){}
		UTMZone(const char& latDiv, const uint8_t& lonDiv):latDiv(latDiv), lonDiv(lonDiv){ }
	};

	struct Ellipse{
		Ellipse(const double& a, const double& f_inv, const char* name) :f_inv(f_inv), a(a), name(name){
		}
		double f_inv;
		double a;
		std::string name;
	};

	class GeoConverter{
	public:
		/*
		! WARN:
			Do not modify the value of this enumerator. 
			The values are associated with a Ellipse instance array.
		*/
		enum EllipseModel{
			WGS84 = 0,
			GRS80 = 1,
			Bessel1841 = 2
		};

		/*
		Input 
			model : the index of EllipseModel. 
		*/
		GeoConverter(EllipseModel model = WGS84){
			selectEllipseModel(model);
		}
		
		/*
		Input 
			model : the index of EllipseModel. 
		*/
		void selectEllipseModel(EllipseModel model);
		
		/*
		Input
			lat : latitude in degrees
			lon : longitude in degrees
		Output
			zone : UTM zone in which (lat, lon) is located
			easting : eastings in meters
			northing : northings in meters
			k : scale factor
		*/
		void convertLL2UTM(double lat, double lon, UTMZone& zone, double& easting, double& northing, double& k);

		/*
		*/
		void convertUTM2LL(){
			//TODO
		}

		/*
		Map a (lat, lon) coordinate into a UTM zone and calculate reference longitude.

		Input
			lat : latitude in degrees
			lon : longitude in degrees
		Output
			refLon : reference longitude value of current longitude zone(1~60). This value is used in UTM coordinates calculation.
			zone : UTM zone in which (lat, lon) is located
		*/
		static bool findZone(double lat,double lon, double& refLon, UTMZone& zone);
		
		/*
		Return currently selected Ellipse instance.

		Return
			a Ellipse instance.
		*/
		const Ellipse& getEllipseModel() const;
	protected:
		double deg2rad(const double& deg){
			return deg * 0.01745329251;
		}
		double rad2deg(const double& rad){
			return rad * 57.29577951308;
		}
	protected:
		EllipseModel ellipse_model;
		double e2, e4, e6, he2;

	};
}
#endif