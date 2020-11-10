
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

Rcpp::List csFs(std::vector<double> uIn, std::vector<double> vIn, std::vector<double> wIn, const double& ZIn, const double& ustarIn, const double& LinvIn, const double& ZoIn, const double& bwIn, const double& sUustarIn, const double& sVustarIn, const double& kvIn, const double& C0In, const double& alphaIn, const double& MaxFetchIn)
	{
	Rcpp::RNGScope scope;

	Rcpp::NumericVector Ri(1003);
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
	double Time, u, v, w, x, y, z, U, dUdz, bsquare, deltaT, deltaXx, deltaXy, deltaXz, fracZ;

	
	const int N = uIn.size();
	int i;
	
	for(int ID = 0; ID < N; ID++){
		x = 0.0;
		y = 0.0;
		z = ZIn;
		u = uIn[ID];
		v = vIn[ID];
		w = wIn[ID];
		Time = 0.0;
		Ri = Rcpp::rnorm(1003);
		i = 0;
		while((z < 1000.0) && (x > -MaxFetchIn)){
			U = ukv*(std::log(z/ZoIn) + z*dpsiMdz - psiMZo);
			dUdz = ukv*(1.0/z + dpsiMdz);
			bsquare = cu3kv*(1.0/z + 5.0*LinvIn);
			deltaT = -alpha2sW2/bsquare;

			u += (s2inv*bsquare*(sigmaW2*(u - U) + ustar2*w) + w*dUdz)*deltaT + Ri[i]*bisqdT;
			i += 1;
			v += bsquare*v*sigmaV22inv*deltaT + Ri[i]*bisqdT;
			i += 1;
			w += s2inv*bsquare*(ustar2*(u - U) + sigmaU2*w)*deltaT + Ri[i]*bisqdT;
			i += 1;
			if(i >= 1000){
				Ri = Rcpp::rnorm(1003);
				i = 0;
			}

			deltaXx = u*deltaT;
			deltaXy = v*deltaT;
			deltaXz = w*deltaT;
			
			if((z+deltaXz)<ZoIn){
				fracZ = (ZoIn - z)/deltaXz;
				x += deltaXx*fracZ;
				y += deltaXy*fracZ;
				Time += deltaT*fracZ;

				IDOut.push_back(ID);
				xOut.push_back(x);
				yOut.push_back(y);
				TimeOut.push_back(Time);
				wTDOut.push_back(w);

				u = 2.0*U - u;
				v = -v;
				w = -w;

				x += u*deltaT*(1.0 - fracZ);
				y -= deltaXy*(1.0 - fracZ);
				z += 2.0*(ZoIn - z) - deltaXz;
				Time += deltaT*(1.0 - fracZ);				
			} else {
				x += deltaXx;
				y += deltaXy;
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

Rcpp::List csFi(std::vector<double> uIn, std::vector<double> vIn, std::vector<double> wIn, const double& ZIn, const double& ustarIn, const double& LinvIn, const double& ZoIn, const double& bwIn, const double& sUustarIn, const double& sVustarIn, const double& kvIn, const double& C0In, const double& alphaIn, const double& MaxFetchIn)
	{
	Rcpp::RNGScope scope;


	Rcpp::NumericVector Ri(1003);
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
	double Time, u, v, w, x, y, z, zL, xi, xi2, psiM, phiE, dsW2dz, sigmaW2, s2inv, U, dUdz, powW, bsquare, deltaT, bisqdT, deltaXx, deltaXy, deltaXz, fracZ;


	const int N = uIn.size();
	int i;


	for(int ID = 0; ID < N; ID++){
		x = 0.0;
		y = 0.0;
		z = ZIn;
		u = uIn[ID];
		v = vIn[ID];
		w = wIn[ID];
		Time = 0.0;
		Ri = Rcpp::rnorm(1003);
		i = 0;

		while((z < 1000.0) && (x > -MaxFetchIn)){
			zL = z*LinvIn;
			xi2 = std::sqrt(1.0-16.0*zL);
			xi = std::sqrt(xi2);
			psiM = std::log(8.0/(1.0+2.0*xi+xi2)/(1.0+xi2)) + 2.0*std::atan(xi) - PI*0.5;
			powW = std::pow(1.0 - 3.0*zL,1.0/3.0);
			phiE = (bw2*bw2*(1.0 - 3.0*zL) + 1.0/powW)/((bw2*bw2 + 1.0)*std::sqrt(std::sqrt(1.0 - 6.0*zL)));
			dsW2dz = -2.0*LinvIn/powW*bw2*ustar2;

			sigmaW2 = ustar2*powW*powW*bw2;
			s2inv = 0.5/(sigmaU2*sigmaW2 - ustar4);

			U = ukv*(std::log(z/ZoIn) + psiM - psiMZo);
			dUdz = ukv/z/xi;

			bsquare = cu3kv*phiE/z;
			deltaT = -alpha2*sigmaW2/bsquare;
			bisqdT = std::sqrt(alpha2*sigmaW2);

			u += (s2inv*bsquare*(sigmaW2*(u - U) + ustar2*w) + w*dUdz)*deltaT + Ri[i]*bisqdT;
			i += 1;
			v += bsquare*v*sigmaV22inv*deltaT + Ri[i]*bisqdT;
			i += 1;
			w += (s2inv*bsquare*(ustar2*(u - U) + sigmaU2*w) + dsW2dz*(0.5 + s2inv*(ustar2*(u - U)*w + sigmaU2*w*w)))*deltaT + Ri[i]*bisqdT;
			i += 1;
			if(i >= 1000){
				Ri = Rcpp::rnorm(1003);
				i = 0;
			}

			deltaXx = u*deltaT;
			deltaXy = v*deltaT;
			deltaXz = w*deltaT;

			if((z+deltaXz)<ZoIn){
				fracZ = (ZoIn - z)/deltaXz;
				x += deltaXx*fracZ;
				y += deltaXy*fracZ;
				Time += deltaT*fracZ;

				IDOut.push_back(ID);
				xOut.push_back(x);
				yOut.push_back(y);
				TimeOut.push_back(Time);
				wTDOut.push_back(w);

				u = 2.0*U - u;
				v = -v;
				w = -w;

				x += u*deltaT*(1.0 - fracZ);
				y -= deltaXy*(1.0 - fracZ);
				z += 2.0*(ZoIn - z) - deltaXz;
				Time += deltaT*(1.0 - fracZ);				
			} else {
				x += deltaXx;
				y += deltaXy;
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

Rcpp::List csFsRK4(std::vector<double> uIn, std::vector<double> vIn, std::vector<double> wIn, const double& ZIn, const double& ustarIn, const double& LinvIn, const double& ZoIn, const double& bwIn, const double& sUustarIn, const double& sVustarIn, const double& kvIn, const double& C0In, const double& alphaIn, const double& MaxFetchIn)
	{
	Rcpp::RNGScope scope;

	Rcpp::NumericVector Ri(1003);
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
	double Time, u, v, w, x, y, z, U, dUdz, bsquare, deltaT, deltaXx, deltaXy, deltaXz, fracZ;
	double h0, k0u, k0v, k0w, l0x, l0y, l0z;
	double z1, u1, v1, w1, h1, k1u, k1v, k1w, l1x, l1y, l1z;
	double z2, u2, v2, w2, h2, k2u, k2v, k2w, l2x, l2y, l2z;
	double z3, u3, v3, w3, h3, k3u, k3v, k3w, l3x, l3y, l3z;

	
	const int N = uIn.size();
	int i;
	
	for(int ID = 0; ID < N; ID++){
		x = 0.0;
		y = 0.0;
		z = ZIn;
		u = uIn[ID];
		v = vIn[ID];
		w = wIn[ID];
		Time = 0.0;
		Ri = Rcpp::rnorm(1003);
		i = 0;
		while((z < 1000.0) && (x > -MaxFetchIn)){
			// k0:
			U = ukv*(std::log(z/ZoIn) + z*dpsiMdz - psiMZo);
			dUdz = ukv*(1.0/z + dpsiMdz);
			bsquare = cu3kv*(1.0/z + 5.0*LinvIn);
			h0 = -alpha2sW2/bsquare;

			k0u = (s2inv*bsquare*(sigmaW2*(u - U) + ustar2*w) + w*dUdz)*h0 + Ri[i]*bisqdT;
			k0v = bsquare*v*sigmaV22inv*h0 + Ri[i+1]*bisqdT;
			k0w = s2inv*bsquare*(ustar2*(u - U) + sigmaU2*w)*h0 + Ri[i+2]*bisqdT;
			
			l0x = u*h0;
			l0y = v*h0;
			l0z = w*h0;

			// k1:
			z1 = ZoIn + std::abs(z + l0z/2.0 - ZoIn); // ev l0z als l0z/2 absp
			u1 = u + k0u/2.0;
			v1 = v + k0v/2.0;
			w1 = w + k0w/2.0;

			U = ukv*(std::log(z1/ZoIn) + z1*dpsiMdz - psiMZo);
			dUdz = ukv*(1.0/z1 + dpsiMdz);
			bsquare = cu3kv*(1.0/z1 + 5.0*LinvIn);

			h1 = -alpha2sW2/bsquare;

			k1u = (s2inv*bsquare*(sigmaW2*(u1 - U) + ustar2*w1) + w1*dUdz)*h1 + Ri[i]*bisqdT;
			k1v = bsquare*v1*sigmaV22inv*h1 + Ri[i+1]*bisqdT;
			k1w = s2inv*bsquare*(ustar2*(u1 - U) + sigmaU2*w1)*h1 + Ri[i+2]*bisqdT;
			
			l1x = u1*h1;
			l1y = v1*h1;
			l1z = w1*h1;

			// k2:
			z2 = ZoIn + std::abs(z + l1z/2.0 - ZoIn);
			u2 = u + k1u/2.0;
			v2 = v + k1v/2.0;
			w2 = w + k1w/2.0;

			U = ukv*(std::log(z2/ZoIn) + z2*dpsiMdz - psiMZo);
			dUdz = ukv*(1.0/z2 + dpsiMdz);
			bsquare = cu3kv*(1.0/z2 + 5.0*LinvIn);
			h2 = -alpha2sW2/bsquare;

			k2u = (s2inv*bsquare*(sigmaW2*(u2 - U) + ustar2*w2) + w2*dUdz)*h2 + Ri[i]*bisqdT;
			k2v = bsquare*v2*sigmaV22inv*h2 + Ri[i+1]*bisqdT;
			k2w = s2inv*bsquare*(ustar2*(u2 - U) + sigmaU2*w2)*h2 + Ri[i+2]*bisqdT;

			l2x = u2*h2;
			l2y = v2*h2;
			l2z = w2*h2;

			// k3:
			z3 = ZoIn + std::abs(z + l2z - ZoIn);
			u3 = u + k2u;
			v3 = v + k2v;
			w3 = w + k2w;

			U = ukv*(std::log(z3/ZoIn) + z3*dpsiMdz - psiMZo);
			dUdz = ukv*(1.0/z3 + dpsiMdz);
			bsquare = cu3kv*(1.0/z3 + 5.0*LinvIn);
			h3 = -alpha2sW2/bsquare;

			k3u = (s2inv*bsquare*(sigmaW2*(u3 - U) + ustar2*w3) + w3*dUdz)*h3 + Ri[i]*bisqdT;
			k3v = bsquare*v3*sigmaV22inv*h3 + Ri[i+1]*bisqdT;
			k3w = s2inv*bsquare*(ustar2*(u3 - U) + sigmaU2*w3)*h3 + Ri[i+2]*bisqdT;
			
			l3x = u3*h3;
			l3y = v3*h3;
			l3z = w3*h3;

			// final steps:

			u += (k0u + (k1u + k2u)*2.0 + k3u)/6.0;
			v += (k0v + (k1v + k2v)*2.0 + k3v)/6.0;
			w += (k0w + (k1w + k2w)*2.0 + k3w)/6.0;
			i += 3;
			if(i >= 1000){
				Ri = Rcpp::rnorm(1003);
				i = 0;
			}

			deltaXx = (l0x + (l1x + l2x)*2.0 + l3x)/6.0;
			deltaXy = (l0y + (l1y + l2y)*2.0 + l3y)/6.0;
			deltaXz = (l0z + (l1z + l2z)*2.0 + l3z)/6.0;

			deltaT = (h0 + (h1 + h2)*2.0 + h3)/6.0;

			if((z+deltaXz)<ZoIn){
				fracZ = (ZoIn - z)/deltaXz;
				x += deltaXx*fracZ;
				y += deltaXy*fracZ;
				Time += deltaT*fracZ;

				IDOut.push_back(ID);
				xOut.push_back(x);
				yOut.push_back(y);
				TimeOut.push_back(Time);
				wTDOut.push_back(deltaXz/deltaT);

				v = -v;
				w = -w;

				z += 2.0*(ZoIn - z) - deltaXz;
				U = ukv*(std::log(z/ZoIn) + z*dpsiMdz - psiMZo);

				u = 2.0*U - u;
				x += deltaT*(1.0 - fracZ)*(2.0*U - deltaXx/deltaT);
				y -= deltaXy*(1.0 - fracZ);
				Time += deltaT*(1.0 - fracZ);				
			} else {
				x += deltaXx;
				y += deltaXy;
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

Rcpp::List csFiRK4(std::vector<double> uIn, std::vector<double> vIn, std::vector<double> wIn, const double& ZIn, const double& ustarIn, const double& LinvIn, const double& ZoIn, const double& bwIn, const double& sUustarIn, const double& sVustarIn, const double& kvIn, const double& C0In, const double& alphaIn, const double& MaxFetchIn)
	{
	Rcpp::RNGScope scope;


	Rcpp::NumericVector Ri(1003);
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
	double Time, u, v, w, x, y, z, zL, xi, xi2, psiM, phiE, dsW2dz, sigmaW2, s2inv, U, dUdz, powW, bsquare, deltaT, bisqdT, deltaXx, deltaXy, deltaXz, fracZ;
	double h0, k0u, k0v, k0w, l0x, l0y, l0z;
	double z1, u1, v1, w1, h1, k1u, k1v, k1w, l1x, l1y, l1z;
	double z2, u2, v2, w2, h2, k2u, k2v, k2w, l2x, l2y, l2z;
	double z3, u3, v3, w3, h3, k3u, k3v, k3w, l3x, l3y, l3z;


	const int N = uIn.size();
	int i;


	for(int ID = 0; ID < N; ID++){
		x = 0.0;
		y = 0.0;
		z = ZIn;
		u = uIn[ID];
		v = vIn[ID];
		w = wIn[ID];
		Time = 0.0;
		Ri = Rcpp::rnorm(1003);
		i = 0;

		while((z < 1000.0) && (x > -MaxFetchIn)){
			// k0:
			zL = z*LinvIn;
			xi2 = std::sqrt(1.0-16.0*zL);
			xi = std::sqrt(xi2);
			psiM = std::log(8.0/(1.0+2.0*xi+xi2)/(1.0+xi2)) + 2.0*std::atan(xi) - PI*0.5;
			powW = std::pow(1.0 - 3.0*zL,1.0/3.0);
			phiE = (bw2*bw2*(1.0 - 3.0*zL) + 1.0/powW)/((bw2*bw2 + 1.0)*std::sqrt(std::sqrt(1.0 - 6.0*zL)));
			dsW2dz = -2.0*LinvIn/powW*bw2*ustar2;

			sigmaW2 = ustar2*powW*powW*bw2;
			s2inv = 0.5/(sigmaU2*sigmaW2 - ustar4);

			U = ukv*(std::log(z/ZoIn) + psiM - psiMZo);
			dUdz = ukv/z/xi;

			bsquare = cu3kv*phiE/z;
			h0 = -alpha2*sigmaW2/bsquare;
			bisqdT = std::sqrt(-bsquare*h0);

			k0u = (s2inv*bsquare*(sigmaW2*(u - U) + ustar2*w) + w*dUdz)*h0 + Ri[i]*bisqdT;
			k0v = bsquare*v*sigmaV22inv*h0 + Ri[i+1]*bisqdT;
			k0w = (s2inv*bsquare*(ustar2*(u - U) + sigmaU2*w) + dsW2dz*(0.5 + s2inv*(ustar2*(u - U)*w + sigmaU2*w*w)))*h0 + Ri[i+2]*bisqdT;
			
			l0x = u*h0;
			l0y = v*h0;
			l0z = w*h0;

			// k1:
			z1 = ZoIn + std::abs(z + l0z/2.0 - ZoIn); // ev l0z als l0z/2 absp
			u1 = u + k0u/2.0;
			v1 = v + k0v/2.0;
			w1 = w + k0w/2.0;
			zL = z1*LinvIn;
			xi2 = std::sqrt(1.0-16.0*zL);
			xi = std::sqrt(xi2);
			psiM = std::log(8.0/(1.0+2.0*xi+xi2)/(1.0+xi2)) + 2.0*std::atan(xi) - PI*0.5;
			powW = std::pow(1.0 - 3.0*zL,1.0/3.0);
			phiE = (bw2*bw2*(1.0 - 3.0*zL) + 1.0/powW)/((bw2*bw2 + 1.0)*std::sqrt(std::sqrt(1.0 - 6.0*zL)));
			dsW2dz = -2.0*LinvIn/powW*bw2*ustar2;

			sigmaW2 = ustar2*powW*powW*bw2;
			s2inv = 0.5/(sigmaU2*sigmaW2 - ustar4);

			U = ukv*(std::log(z1/ZoIn) + psiM - psiMZo);
			dUdz = ukv/z1/xi;

			bsquare = cu3kv*phiE/z1;

			h1 = -alpha2*sigmaW2/bsquare;
			bisqdT = std::sqrt(-bsquare*h1);

			k1u = (s2inv*bsquare*(sigmaW2*(u1 - U) + ustar2*w1) + w1*dUdz)*h1 + Ri[i]*bisqdT;
			k1v = bsquare*v1*sigmaV22inv*h1 + Ri[i+1]*bisqdT;
			k1w = (s2inv*bsquare*(ustar2*(u1 - U) + sigmaU2*w1) + dsW2dz*(0.5 + s2inv*(ustar2*(u1 - U)*w1 + sigmaU2*w1*w1)))*h1 + Ri[i+2]*bisqdT;
			
			l1x = u1*h1;
			l1y = v1*h1;
			l1z = w1*h1;

			// k2:
			z2 = ZoIn + std::abs(z + l1z/2.0 - ZoIn);
			u2 = u + k1u/2.0;
			v2 = v + k1v/2.0;
			w2 = w + k1w/2.0;
			zL = z2*LinvIn;
			xi2 = std::sqrt(1.0-16.0*zL);
			xi = std::sqrt(xi2);
			psiM = std::log(8.0/(1.0+2.0*xi+xi2)/(1.0+xi2)) + 2.0*std::atan(xi) - PI*0.5;
			powW = std::pow(1.0 - 3.0*zL,1.0/3.0);
			phiE = (bw2*bw2*(1.0 - 3.0*zL) + 1.0/powW)/((bw2*bw2 + 1.0)*std::sqrt(std::sqrt(1.0 - 6.0*zL)));
			dsW2dz = -2.0*LinvIn/powW*bw2*ustar2;

			sigmaW2 = ustar2*powW*powW*bw2;
			s2inv = 0.5/(sigmaU2*sigmaW2 - ustar4);

			U = ukv*(std::log(z2/ZoIn) + psiM - psiMZo);
			dUdz = ukv/z2/xi;

			bsquare = cu3kv*phiE/z2;
			h2 = -alpha2*sigmaW2/bsquare;
			bisqdT = std::sqrt(-bsquare*h2);

			k2u = (s2inv*bsquare*(sigmaW2*(u2 - U) + ustar2*w2) + w2*dUdz)*h2 + Ri[i]*bisqdT;
			k2v = bsquare*v2*sigmaV22inv*h2 + Ri[i+1]*bisqdT;
			k2w = (s2inv*bsquare*(ustar2*(u2 - U) + sigmaU2*w2) + dsW2dz*(0.5 + s2inv*(ustar2*(u2 - U)*w2 + sigmaU2*w2*w2)))*h2 + Ri[i+2]*bisqdT;

			l2x = u2*h2;
			l2y = v2*h2;
			l2z = w2*h2;

			// k3:
			z3 = ZoIn + std::abs(z + l2z - ZoIn);
			u3 = u + k2u;
			v3 = v + k2v;
			w3 = w + k2w;
			zL = z3*LinvIn;
			xi2 = std::sqrt(1.0-16.0*zL);
			xi = std::sqrt(xi2);
			psiM = std::log(8.0/(1.0+2.0*xi+xi2)/(1.0+xi2)) + 2.0*std::atan(xi) - PI*0.5;
			powW = std::pow(1.0 - 3.0*zL,1.0/3.0);
			phiE = (bw2*bw2*(1.0 - 3.0*zL) + 1.0/powW)/((bw2*bw2 + 1.0)*std::sqrt(std::sqrt(1.0 - 6.0*zL)));
			dsW2dz = -2.0*LinvIn/powW*bw2*ustar2;

			sigmaW2 = ustar2*powW*powW*bw2;
			s2inv = 0.5/(sigmaU2*sigmaW2 - ustar4);

			U = ukv*(std::log(z3/ZoIn) + psiM - psiMZo);
			dUdz = ukv/z3/xi;

			bsquare = cu3kv*phiE/z3;
			h3 = -alpha2*sigmaW2/bsquare;
			bisqdT = std::sqrt(-bsquare*h3);

			k3u = (s2inv*bsquare*(sigmaW2*(u3 - U) + ustar2*w3) + w3*dUdz)*h3 + Ri[i]*bisqdT;
			k3v = bsquare*v3*sigmaV22inv*h3 + Ri[i+1]*bisqdT;
			k3w = (s2inv*bsquare*(ustar2*(u3 - U) + sigmaU2*w3) + dsW2dz*(0.5 + s2inv*(ustar2*(u3 - U)*w3 + sigmaU2*w3*w3)))*h3 + Ri[i+2]*bisqdT;
			
			l3x = u3*h3;
			l3y = v3*h3;
			l3z = w3*h3;

			// final steps:

			u += (k0u + (k1u + k2u)*2.0 + k3u)/6.0;
			v += (k0v + (k1v + k2v)*2.0 + k3v)/6.0;
			w += (k0w + (k1w + k2w)*2.0 + k3w)/6.0;
			i += 3;
			if(i >= 1000){
				Ri = Rcpp::rnorm(1003);
				i = 0;
			}

			deltaXx = (l0x + (l1x + l2x)*2.0 + l3x)/6.0;
			deltaXy = (l0y + (l1y + l2y)*2.0 + l3y)/6.0;
			deltaXz = (l0z + (l1z + l2z)*2.0 + l3z)/6.0;

			deltaT = (h0 + (h1 + h2)*2.0 + h3)/6.0;

			if((z+deltaXz)<ZoIn){
				fracZ = (ZoIn - z)/deltaXz;
				x += deltaXx*fracZ;
				y += deltaXy*fracZ;
				Time += deltaT*fracZ;

				IDOut.push_back(ID);
				xOut.push_back(x);
				yOut.push_back(y);
				TimeOut.push_back(Time);
				wTDOut.push_back(deltaXz/deltaT);

				v = -v;
				w = -w;

				z += 2.0*(ZoIn - z) - deltaXz;
				zL = z*LinvIn;
				xi2 = std::sqrt(1.0-16.0*zL);
				xi = std::sqrt(xi2);
				psiM = std::log(8.0/(1.0+2.0*xi+xi2)/(1.0+xi2)) + 2.0*std::atan(xi) - PI*0.5;
				U = ukv*(std::log(z/ZoIn) + psiM - psiMZo);

				u = 2.0*U - u;
				x += deltaT*(1.0 - fracZ)*(2.0*U - deltaXx/deltaT);
				y -= deltaXy*(1.0 - fracZ);
				Time += deltaT*(1.0 - fracZ);				
			} else {
				x += deltaXx;
				y += deltaXy;
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

