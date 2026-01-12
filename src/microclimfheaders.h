// microclimfheaders.h
// Radiation model
#include <Rcpp.h>
using namespace Rcpp;
struct solmodel {
    double zend;
    double zenr;
    double azid;
    double azir;
};
struct kstruct {
    double k;
    double kd;
    double Kc;
};
struct tsdifstruct {
    double p1;
    double p2;
    double p3;
    double p4;
    double om;
    double a;
    double gma;
    double J;
    double del;
    double h;
    double u1;
    double S1;
    double D1;
    double D2;
};
struct tsdirstruct {
    double sig;
    double p5;
    double p6;
    double p7;
    double p8;
    double p9;
    double p10;
};
struct radmodel {
    std::vector<double> ground; // ground absorbed radiation (W/m^2)
    std::vector<double> canopy; // canopy absorbed radiation (W/m^2)
    std::vector<double> albedo; // albedo
};
struct tirstruct {
    double albd;
    double Rddn_g;
    double Rdup_z;
    double Rddn_z;
    double gi;
    double trdn;
    double trdu;
    double amx;
    double pait;
    double paiaa;
    double om;
    double omp;
    double a;
    double gma;
    double J;
    double del;
    double h;
    double u1;
    double S1;
    double D1;
    double D2;
};
struct radmodel2 {
    double radGsw;
    double radGlw;
    double radCsw;
    double radClw;
    double Rbdown;
    double Rddown;
    double Rdup;
    double radLsw;
    double radLpar;
    double lwout;
    double zend;
};
struct stompstruct {
    double Rsmx;
    double psiw0;
    double kk;
    double rat;
};
struct Gmodel {
    std::vector<double> G;
    std::vector<double> Gmin;
    std::vector<double> Gmax;
};
struct tiwstruct {
    double d;
    double zm;
    double a;
};
struct windmodel {
    double uf;
    double uz;
    double gHa;
};
struct penmonstruct {
    double Ts;
    double H;
    double L;
    double Rem;
    double mu;
};
struct soilstruct {
    double c1;
    double c3;
    double c4;
};
struct soilkstruct {
    double k;
    double DD;
};
struct soilpstruct {
    double Smax;
    double Smin;
    double soilb;
    double psi_e;
    double Vq;
    double Vm;
    double Mc;
    double rho;
};
struct soilmodelG0 {
    double Tg;
    double Rnet;
    double surfwet;
    double radabs;
};
struct soilmodel {
    double Tg;
    double G;
    double DD;
};
struct abovecanstruct {
    double Tz;
    double ez;
};
struct leaftempstruct {
    double tleaf;
    double H;
    double L;
    double lwdn;
    double lwup;
};
struct abovemodel {
    double Tz;
    double tleaf;
    double rh;
    double lwdn;
    double lwup;
};
// Used in snow model
struct obspoint {
    int year;
    int month;
    int day;
    double hour;
};
struct climpoint {
    double tc;
    double ea;
    double pk;
    double u2;
    double Rsw;
    double Rdif;
    double Rlw;
    double prec;
    double Tci; // temperature of snow pack surface
    double te; // estimated snow surface temperature
};
struct vegpoint {
    double pai;
    double hgt;
    double ltra;
    double clump;
};
struct snowpoint {
    double sdenc;
    double sdeng;
    double sdepc;
    double sdepg;
    double snowagec;
    double snowageg;
    double alb;
};
struct snowpoint2 {
    double snowtempg; 
    double snowtempc;
    double sdepc;
    double sdepg;
    double sdenc;
    double albg;
    double albc;
};
struct otherpoint {
    double slope;
    double aspect;
    double lat;
    double lon;
    double zref;
    double psim;
    double psih;
    double G;
};
struct snowrad {
    double RabsC;
    double RswabsG;
    double RlwabsG;
    double tr; // transmission through canopy
};
struct snowmodpoint {
    // Canopy + ground
    double Tc;
    double mSc; // sublimation (m SWE)
    double mMc; // temperature melt (m SWE)
    double mRc; // rain melt (m SWE)
    // Ground
    double Tg;
    double mSg; // sublimation (m SWE)
    double mMg; // temperature melt (m SWE)
    double mRg; // rain melt (m SWE)
    // Other
    double cis; // canopy interception
    double uf; // fraction velocity
    double RswabsG; // shortwave radiation absorbed by ground
    double RlwabsG; // longwave radiation absorbed by ground
    double tr; // tranmission through canopy
    double gHa; // Convective conductance
    double Tcp; // snow surface temperature (allowing < 0 deg c)
    // New parameters inherted form snow
    double snowagec; // age of canopy snow layer
    double snowageg; // age of ground snow layer
    double sdenc; // density of canopy layer
    double sdeng; // density of ground layer
    double sdepc; // depth of canopy layer
    double sdepg; // depth of ground layer
    double pai; // plant area index adjusted for snow
    double hgt; // vegetation height adjusted for snow
};
struct snowmicro {
    double Tz;
    double tleaf;
    double rh;
    double uz;
    double Rbdown;
    double Rddown;
    double Rlwdn;
    double Rdup;
    double Rlwup;
};

