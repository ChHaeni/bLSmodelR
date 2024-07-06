
#include <Rcpp.h>
using namespace Rcpp;

#ifndef PI 
#define PI 3.14159265
#endif

// [[Rcpp::export]]

Rcpp::List csFs(
        Rcpp::NumericVector uIn, 
        Rcpp::NumericVector vIn, 
        Rcpp::NumericVector wIn, 
        const double& ZIn, 
        const double& ustarIn, 
        const double& LinvIn, 
        const double& ZoIn, 
        const double& bwIn, 
        const double& sUustarIn, 
        const double& sVustarIn, 
        const double& kvIn, 
        const double& C0In, 
        const double& alphaIn, 
        const double& MaxFetchIn
        )
	{
	Rcpp::RNGScope scope;

	/* Rcpp::NumericVector Ri(1003); */
    double Ri;
	std::vector<double> xOut;
	std::vector<double> yOut;
	std::vector<double> wTDOut;
	std::vector<double> TimeOut;
	std::vector<int> IDOut;

	const double ustar2 = ustarIn*ustarIn;
	const double ustar4 = ustar2*ustar2;
	const double sigmaU2 = sUustarIn*sUustarIn*ustar2;
	const double sigmaV22inv = 0.5/(sVustarIn*sVustarIn*ustar2);
	const double bw2 = bwIn*bwIn;
	const double ukv = ustarIn/kvIn;
	const double cu3kv = C0In*ustar2*ukv;
	const double psiMZo = 4.8*ZoIn*LinvIn;
	const double dpsiMdz = 4.8*LinvIn;
	const double sigmaW2 = ustar2*bw2;
	const double s2inv = 0.5/(sigmaU2*sigmaW2 - ustar4);
	const double alpha2sW2 = alphaIn*2.0*sigmaW2;
	const double bisqdT = std::sqrt(alpha2sW2);
	double Time, u, v, w, x, y, z, U, bsquare, deltaT, deltaXz, fracZ;

	
	const int N = uIn.size();
	
	for(int ID = 0; ID < N; ID++){
		x = 0.0;
		y = 0.0;
		z = ZIn;
		u = uIn[ID];
		v = vIn[ID];
		w = wIn[ID];
		Time = 0.0;
		while((z < 1000.0) && (x > -MaxFetchIn)){
            Ri = R::rnorm(0.0, 1.0);
			U = ukv*(std::log(z/ZoIn) + z*dpsiMdz - psiMZo);
			bsquare = cu3kv*(1.0/z + 5.0*LinvIn);
			deltaT = -alpha2sW2/bsquare;

			u += (s2inv*bsquare*(sigmaW2*(u - U) + ustar2*w) + w*ukv*(1.0/z + dpsiMdz))*deltaT + Ri*bisqdT;
			v += bsquare*v*sigmaV22inv*deltaT + Ri*bisqdT;
			w += s2inv*bsquare*(ustar2*(u - U) + sigmaU2*w)*deltaT + Ri*bisqdT;

			deltaXz = w*deltaT;
			
			if((z+deltaXz)<ZoIn){
				fracZ = (ZoIn - z)/deltaXz;
				x += u*deltaT*fracZ;
				y += v*deltaT*fracZ;
				Time += deltaT*fracZ;

				IDOut.push_back(ID + 1);
				xOut.push_back(x);
				yOut.push_back(y);
				TimeOut.push_back(Time);
				wTDOut.push_back(w);

				u = 2.0*U - u;
				v = -v;
				w = -w;

				x += u*deltaT*(1.0 - fracZ);
				y += v*deltaT*(1.0 - fracZ);
				z += 2.0*(ZoIn - z) - deltaXz;
				Time += deltaT*(1.0 - fracZ);				
			} else {
				x += u*deltaT;
				y += v*deltaT;
				z += deltaXz;
				Time += deltaT;				
			}
		}
	}

	return Rcpp::List::create(
		_["Traj_IDOut"] = IDOut,
		_["TimeOut"] = TimeOut,
		_["xOut"] = xOut,
		_["yOut"] = yOut,
		_["wTDOut"] = wTDOut);
}


// [[Rcpp::export]]

Rcpp::List csFi(
        Rcpp::NumericVector uIn, 
        Rcpp::NumericVector vIn, 
        Rcpp::NumericVector wIn, 
        const double& ZIn, 
        const double& ustarIn, 
        const double& LinvIn, 
        const double& ZoIn, 
        const double& bwIn, 
        const double& sUustarIn, 
        const double& sVustarIn, 
        const double& kvIn, 
        const double& C0In, 
        const double& alphaIn, 
        const double& MaxFetchIn
        )
	{
	Rcpp::RNGScope scope;


	/* Rcpp::NumericVector Ri(1003); */
    double Ri;
	std::vector<double> xOut;
	std::vector<double> yOut;
	std::vector<double> wTDOut;
	std::vector<double> TimeOut;
	std::vector<int> IDOut;


	const double ustar2 = ustarIn*ustarIn;
	const double ustar4 = ustar2*ustar2;
	const double sigmaU2 = sUustarIn*sUustarIn*ustar2;
	const double sigmaV22inv = 0.5/(sVustarIn*sVustarIn*ustar2);
	const double bw2 = bwIn*bwIn;
	const double x2Zo = std::sqrt(1.0-16.0*ZoIn*LinvIn);
	const double xZo = std::sqrt(x2Zo);
	const double ukv = ustarIn/kvIn;
	const double cu3kv = C0In*ustar2*ukv;
	const double psiMZo = std::log(8.0/(1.0+2.0*xZo+x2Zo)/(1.0+x2Zo)) + 2.0*std::atan(xZo) - PI*0.5;
	const double alpha2 = 2.0*alphaIn;
	double Time, u, v, w, x, y, z, zL, xi, xi2, sigmaW2, s2inv, U, powW, bsquare, deltaT, bisqdT, deltaXz, fracZ;


	const int N = uIn.size();


	for(int ID = 0; ID < N; ID++){
		x = 0.0;
		y = 0.0;
		z = ZIn;
		u = uIn[ID];
		v = vIn[ID];
		w = wIn[ID];
		Time = 0.0;

		while((z < 1000.0) && (x > -MaxFetchIn)){
            Ri = R::rnorm(0.0, 1.0);
			zL = z*LinvIn;
			xi2 = std::sqrt(1.0-16.0*zL);
			xi = std::sqrt(xi2);
			powW = std::pow(1.0 - 3.0*zL,1.0/3.0);

			sigmaW2 = ustar2*powW*powW*bw2;
			s2inv = 0.5/(sigmaU2*sigmaW2 - ustar4);

			U = ukv*(std::log(z/ZoIn) + std::log(8.0/(1.0+2.0*xi+xi2)/(1.0+xi2)) + 2.0*std::atan(xi) - PI*0.5 - psiMZo);

			bsquare = cu3kv*(bw2*bw2*(1.0 - 3.0*zL) + 1.0/powW)/((bw2*bw2 + 1.0)*std::sqrt(std::sqrt(1.0 - 6.0*zL)))/z;
			deltaT = -alpha2*sigmaW2/bsquare;
			bisqdT = std::sqrt(alpha2*sigmaW2);

			u += (s2inv*bsquare*(sigmaW2*(u - U) + ustar2*w) + w*ukv/z/xi)*deltaT + Ri*bisqdT;
			v += bsquare*v*sigmaV22inv*deltaT + Ri*bisqdT;
			w += (s2inv*bsquare*(ustar2*(u - U) + sigmaU2*w) - 2.0*LinvIn/powW*bw2*ustar2*(0.5 + s2inv*(ustar2*(u - U)*w + sigmaU2*w*w)))*deltaT + Ri*bisqdT;

			deltaXz = w*deltaT;

			if((z+deltaXz)<ZoIn){
				fracZ = (ZoIn - z)/deltaXz;
				x += u*deltaT*fracZ;
				y += v*deltaT*fracZ;
				Time += deltaT*fracZ;

				IDOut.push_back(ID + 1);
				xOut.push_back(x);
				yOut.push_back(y);
				TimeOut.push_back(Time);
				wTDOut.push_back(w);

				u = 2.0*U - u;
				v = -v;
				w = -w;

				x += u*deltaT*(1.0 - fracZ);
				y += v*deltaT*(1.0 - fracZ);
				z += 2.0*(ZoIn - z) - deltaXz;
				Time += deltaT*(1.0 - fracZ);				
			} else {
				x += u*deltaT;
				y += v*deltaT;
				z += deltaXz;
				Time += deltaT;				
			}
		}
	}

	return Rcpp::List::create(
		_["Traj_IDOut"] = IDOut,
		_["TimeOut"] = TimeOut,
		_["xOut"] = xOut,
		_["yOut"] = yOut,
		_["wTDOut"] = wTDOut);
}

