/*================================================================================
  ================================================================================
    ESTE PRGRAMA INYECTA UNA NUBE ESFERICA DE  ELECTRONES 
    EN UN DISPOSITIVO DE GYRAC.

    IMPRIME LA DISTRIBUCIÓN ESPACIAL CON EL CORRESPONDIENTE 
    VALOR DE GAMMA, EN DIFERENTES  (5).
  ================================================================================
  ================================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <iostream>


#define Ne 1000      //numero de electrones           
#define me 9.10938356e-31 // Rest Mass electron (Kg)
#define q -1.6021766208e-19  // Electric charge particle(C)(Electron)
#define c  299792458    // Light Velocity  (m/s)
#define p 1         // Index of TE11p mode
#define npz 100
#define npr 65
#define m 250       // Number Step by Hf Cycle
#define kt 750*310*5     // Total number Step
#define N 3     //numero de campos electricos leidos
#define eps0 8.8541878176e-12  // Electric permitivity
#define npe 1000	//numero de puntos espectro energetico

using namespace std;

//***********************************************
//     FUNCIONES 
//***********************************************

void READ_FILES(void);
void PRINT_FIELD(void);
void PRINT_POSITION(void);
void MOTION(void);
void PHASE_SHIFT(void);
void PHASE_VELOCITY(void);
void PHASE_ELECTRIC_FIELD(void);
void LINEAR_INTERPOLATION(void);
void INJECT_ELECTRON(void);
double distrib_vel(void);

//***********************************************
//      DECLARACIÓN DE VARIABLES 
//***********************************************


double Uz0,shr,shz;
double om,B0,rl,W0,gama,vz0,g0,T,dt,dt2,ang,ang1,ang2;
double gx[m],gy[m],Bx_hf[m],By_hf[m],x[Ne],y[Ne],z[Ne],Ux[Ne],Uy[Ne],Uz[Ne],gmn,gm[Ne];
//double x0=1.0; //cm
double ang1_new,ang1_old,dif1,n1,ang2_new,ang2_old,dif2,n2;
double Ex,Ey;
float Br[npr][npz],Bz[npr][npz];
float Brr[N][npr][npz],Bzr[N][npr][npz];
int kp,kk;
double tm;
double j0( double xx);
double j1 (double xx);
double S11=1.841,h,J1,dJ1;
double sum_energy;
double Bxz_2;
double E0=1,rc=4.54,Lc=10.0; //valores ingresados para este problema
double f=2.45,aK0=1;                   //valores ingresados para este problema [GHz]
double r_c,aL;  //variables de normalizacion 
double x0=0.0001; //cm
double y00=0.001;
double z0=(Lc/2);
double ri,zi;  // Posicion del ion en (r,z) para usar en la interpolacion.
double brp,bzp,bxp,byp; // componenete del campo electrico 
double ene_cont[npe];	//contador de electrones dependiendo de la energía
int pr,pz;  // Indices de malla para la posicion del ion.
int k=0,ik=0;  
int check;  //conteo de electrones inyectados
int cont_e_gm=0,e_gm[Ne];   //conteo de electrones que alcanza alta energía
int cont_e_out=0,e_out[Ne];   //conteo de electrones que salen de la cavidad
int gm_max=5;  //gamma "maximo" el cual sirve de referencia para conteo  
int ie;     //variable para recorer los electrones
int ref;	//referencia para ciclo que decrece en el tiempo
double gmax=1,g_max[npe];	//valores para distribución energetica


//==================== velocidad de inyección ===========
long double Te0=2.0; // (Ti-->kB*Ti [eV]);
long double Te=2;//100.*1000*1.6022e-19; // (Temperature Final (estimada) in KeV
double vte=sqrt(Te0*1.6e-19/me); // electron thermal velocity (m/s) ;
double lambda_Debye=sqrt((eps0*Te)/(1e6*Ne*pow(-q,2)));

//==================== crear muchos archivos de salida ===========
char buffer[80]={0}; // lo usaremos para guardar el nombre del fichero


FILE *out2,*out3,*out4,*out5,*out6,*out7,*out8,*out11,*out1;
FILE *salida1,*salida2,*salida3,*salida4,*salida5,*salida6,*salida7,*salida8,*salida9,*salida10;
FILE *salida11,*salida12,*salida13,*salida14,*salida15,*salida16,*salida17,*salida18,*salida19,*salida20;
FILE *salida21,*salida22,*salida23,*salida24,*salida25,*salida26,*salida27,*salida28,*salida29,*salida30;
FILE *salida31,*salida32,*salida33,*salida34,*salida35,*salida36,*salida37,*salida38,*salida39,*salida40;
FILE *salida41,*salida42,*salida43,*salida44,*salida45,*salida46,*salida47,*salida48,*salida49,*salida50;
FILE *salida51,*salida52,*salida53,*salida54,*salida55,*salida56,*salida57,*salida58,*salida59,*salida60;
FILE *salida61,*salida62,*salida63,*salida64,*salida65,*salida66,*salida67,*salida68,*salida69,*salida70;
FILE *salida71,*salida72,*salida73,*salida74,*salida75,*salida76,*salida77,*salida78,*salida79,*salida80;
FILE *salida81,*salida82,*salida83,*salida84,*salida85,*salida86,*salida87,*salida88,*salida89,*salida90;
FILE *salida91,*salida92,*salida93,*salida94,*salida95,*salida96,*salida97,*salida98,*salida99,*salida100;

int main()
{
    out2=fopen("B_rz_i.txt","w");
    out3=fopen("B_rz_f.txt","w");
    out4=fopen("Posiciones_0.txt","w");
    out5=fopen("Posiciones_25.txt","w");
    out6=fopen("Posiciones_50.txt","w");
    out7=fopen("Posiciones_75.txt","w");
    out8=fopen("Posiciones_100.txt","w");
    out11=fopen("B_z_i_f.txt","w");

    salida1=fopen("resultados1.txt","w");
    salida2=fopen("resultados2.txt","w");
    salida3=fopen("resultados3.txt","w");
    salida4=fopen("resultados4.txt","w");
    salida5=fopen("resultados5.txt","w");
    salida6=fopen("resultados6.txt","w");
    salida7=fopen("resultados7.txt","w");
    salida8=fopen("resultados8.txt","w");
    salida9=fopen("resultados9.txt","w");
    salida10=fopen("resultados10.txt","w");
    salida11=fopen("resultados11.txt","w");
    salida12=fopen("resultados12.txt","w");
    salida13=fopen("resultados13.txt","w");
    salida14=fopen("resultados14.txt","w");
    salida15=fopen("resultados15.txt","w");
    salida16=fopen("resultados16.txt","w");
    salida17=fopen("resultados17.txt","w");
    salida18=fopen("resultados18.txt","w");
    salida19=fopen("resultados19.txt","w");
    salida20=fopen("resultados20.txt","w");   
    salida21=fopen("resultados21.txt","w");
    salida22=fopen("resultados22.txt","w");
    salida23=fopen("resultados23.txt","w");
    salida24=fopen("resultados24.txt","w");
    salida25=fopen("resultados25.txt","w");
    salida26=fopen("resultados26.txt","w");
    salida27=fopen("resultados27.txt","w");
    salida28=fopen("resultados28.txt","w");
    salida29=fopen("resultados29.txt","w");     
    salida30=fopen("resultados30.txt","w");
    salida31=fopen("resultados31.txt","w");
    salida32=fopen("resultados32.txt","w");
    salida33=fopen("resultados33.txt","w");
    salida34=fopen("resultados34.txt","w");
    salida35=fopen("resultados35.txt","w");
    salida36=fopen("resultados36.txt","w");
    salida37=fopen("resultados37.txt","w");
    salida38=fopen("resultados38.txt","w");
    salida39=fopen("resultados39.txt","w");
    salida40=fopen("resultados40.txt","w");    
    salida41=fopen("resultados41.txt","w");
    salida42=fopen("resultados42.txt","w");
    salida43=fopen("resultados43.txt","w");
    salida44=fopen("resultados44.txt","w");
    salida45=fopen("resultados45.txt","w");
    salida46=fopen("resultados46.txt","w");
    salida47=fopen("resultados47.txt","w");
    salida48=fopen("resultados48.txt","w");
    salida49=fopen("resultados49.txt","w");
    salida50=fopen("resultados50.txt","w");
    salida51=fopen("resultados51.txt","w");
    salida52=fopen("resultados52.txt","w");
    salida53=fopen("resultados53.txt","w");
    salida54=fopen("resultados54.txt","w");
    salida55=fopen("resultados55.txt","w");
    salida56=fopen("resultados56.txt","w");
    salida57=fopen("resultados57.txt","w");
    salida58=fopen("resultados58.txt","w");
    salida59=fopen("resultados59.txt","w");
    salida60=fopen("resultados60.txt","w");
    salida61=fopen("resultados61.txt","w");
    salida62=fopen("resultados62.txt","w");
    salida63=fopen("resultados63.txt","w");
    salida64=fopen("resultados64.txt","w");
    salida65=fopen("resultados65.txt","w");
    salida66=fopen("resultados66.txt","w");
    salida67=fopen("resultados67.txt","w");
    salida68=fopen("resultados68.txt","w");
    salida69=fopen("resultados69.txt","w");
    salida70=fopen("resultados70.txt","w");
    salida71=fopen("resultados71.txt","w");
    salida72=fopen("resultados72.txt","w");
    salida73=fopen("resultados73.txt","w");
    salida74=fopen("resultados74.txt","w");
    salida75=fopen("resultados75.txt","w");
    salida76=fopen("resultados76.txt","w");
    salida77=fopen("resultados77.txt","w");
    salida78=fopen("resultados78.txt","w");
    salida79=fopen("resultados79.txt","w");
    salida80=fopen("resultados80.txt","w");
    salida81=fopen("resultados81.txt","w");
    salida82=fopen("resultados82.txt","w");
    salida83=fopen("resultados83.txt","w");
    salida84=fopen("resultados84.txt","w");
    salida85=fopen("resultados85.txt","w");
    salida86=fopen("resultados86.txt","w");
    salida87=fopen("resultados87.txt","w");
    salida88=fopen("resultados88.txt","w");
    salida89=fopen("resultados89.txt","w");
    salida90=fopen("resultados90.txt","w");
    salida91=fopen("resultados91.txt","w");
    salida92=fopen("resultados92.txt","w");
    salida93=fopen("resultados93.txt","w");
    salida94=fopen("resultados94.txt","w");
    salida95=fopen("resultados95.txt","w");
    salida96=fopen("resultados96.txt","w");
    salida97=fopen("resultados97.txt","w");
    salida98=fopen("resultados98.txt","w");
    salida99=fopen("resultados99.txt","w");
    salida100=fopen("Electron_test.txt","w");

    //***********************************************
//      CONVERSIÓN DE UNIDADES AL S.I
//***********************************************


    f=f*1e9;                 //  Frequency (Hz)
    om=2.*M_PI*f;
    B0=fabs(me*om/q);          //  Magnetic field value of clasic resonance.
    rl=c/om;                   //  Relativistic Larmor Radius
    W0=me*pow(c,2);            //  Rest Energy
    aK0=aK0*1.e3*1.6e-19;      //  Energy injection (J)
    gama=1.0+aK0/W0;           //  Initial Relativistic factor
    rc=rc/100.;                //  Radii cavity (m)
    Lc=Lc/100.;                //  Length cavity (m)
    E0=E0*1.e5;                //  (V/m)
    kp=0;
    n1=0;
    n2=0;
    ang1_old=0;
    ang2_old=0;
    x0=x0/100.;  
    y00=y00/100.;
    z0=z0/100.;

    shr=rc/(npr-1);
    shz=Lc/(npz-1);

//***********************************************
//      NPRMALIZACIÓN DE LAS VARIABLES
//***********************************************

    r_c=rc/rl;
    aL=Lc/rl;
    h=S11/r_c;
    g0=E0/(-B0*c);
    vz0=sqrt(1.-pow(gama,-2));  // Injection velocity (in c units)
    Uz0=gama*vz0;               // Initial momentum
    T=2.0*M_PI;
    dt=T/m;
    dt2=0.5*dt;
    x0=x0/rl;
    y00=y00/rl;
    z0=z0/rl;

//***********************************************
//      IMPRIME CONDICIONES INICIALES RELEVANTES
//***********************************************
    printf("Condiciones del sistema \n");

    printf("    E=%e \n", E0);
    
    printf("    om=%e\n",om);

    printf("    B0=%e\n",B0);

    printf("    g0=%e\n",-g0);


    INJECT_ELECTRON();  printf("\nSe inyectaron %d electrones\n",check ); // Inject Plasma.

    READ_FILES();	//lectura de todos los campo

    PRINT_FIELD();  //verifica los campos leidos 


//***********************************************
//      NORMALIZACION DEL CAMPO MAGNETICO -B0
//***********************************************

    for(int ii=1;ii<N;ii++) 
        for(int i=0;i<npr;i++)   
		   for(int j=0;j<npz;j++){
    		     Brr[ii][i][j]=-Brr[ii][i][j];
    		     Bzr[ii][i][j]=-Bzr[ii][i][j];
            }

//***********************************************
//      INICIALIZACION DE ETIQUETAS
//***********************************************

    for (int i = 0; i < Ne; ++i){  
        e_gm[i]=0; //inicializacion para etiquetar particula con gamma "alto"
        e_out[i]=0; //inicializacion para etiquetar particula que sale 
	if(i<=npe){
		ene_cont[i]=0; //inicializacion para conteo energetico
		g_max[i]=gmax;

		}
    }
	      
//***********************************************
//      COMIENZA EL CICLO DEL MOV.
//***********************************************
printf("porcentade de ejecución\n");
for(ik=0;ik<kt;ik++){ // Time Cycle for integration of motion equation

    tm=ik*dt;          // current time
    kk=ik-kp;          // Index for electric field


    for (ie = 0; ie < Ne; ++ie)  MOTION();  //ciclo para que en ik pase por todos los electrones        

    PRINT_POSITION();   //imprime archivos para las posiciones en diferentes momentos 

    if(ik%12000==0)printf(" %f %c \n", ik*pow(kt,-1)*100,'%');

    if (ik==kt-1){
        printf("%d electrones tienen gamma mayor que %d \n", cont_e_gm,gm_max);
        printf("%d electrones salieron de la cavidad \n", cont_e_out);
    }

    if(ik==kp+m-1)
    {
        kp=kp+m;
    }
}

fclose(salida1);
fclose(salida2);
fclose(salida3);
fclose(salida4);
fclose(salida5);
fclose(salida6);
fclose(salida7);
fclose(salida8);
fclose(salida9);
fclose(salida10);
fclose(salida11);
fclose(salida12);
fclose(salida13);
fclose(salida14);
fclose(salida15);
fclose(salida16);
fclose(salida17);
fclose(salida18);
fclose(salida19);
fclose(salida20);
fclose(salida21);
fclose(salida22);
fclose(salida23);
fclose(salida24);
fclose(salida25);
fclose(salida26);
fclose(salida27);
fclose(salida28);
fclose(salida29);
fclose(salida30);
fclose(salida31);
fclose(salida32);
fclose(salida33);
fclose(salida34);
fclose(salida35);
fclose(salida36);
fclose(salida37);
fclose(salida38);
fclose(salida39);
fclose(salida40);
fclose(salida41);
fclose(salida42);
fclose(salida43);
fclose(salida44);
fclose(salida45);
fclose(salida46);
fclose(salida47);
fclose(salida48);
fclose(salida49);
fclose(salida50);
fclose(salida51);
fclose(salida52);
fclose(salida53);
fclose(salida54);
fclose(salida55);
fclose(salida56);
fclose(salida57);
fclose(salida58);
fclose(salida59);
fclose(salida60);
fclose(salida61);
fclose(salida62);
fclose(salida63);
fclose(salida64);
fclose(salida65);
fclose(salida66);
fclose(salida67);
fclose(salida68);
fclose(salida69);
fclose(salida70);
fclose(salida71);
fclose(salida72);
fclose(salida73);
fclose(salida74);
fclose(salida75);
fclose(salida76);
fclose(salida77);
fclose(salida78);
fclose(salida79);
fclose(salida80);
fclose(salida81);
fclose(salida82);
fclose(salida83);
fclose(salida84);
fclose(salida85);
fclose(salida86);
fclose(salida87);
fclose(salida88);
fclose(salida89);
fclose(salida90);
fclose(salida91);
fclose(salida92);
fclose(salida93);
fclose(salida94);
fclose(salida95);
fclose(salida96);
fclose(salida97);
fclose(salida98);
fclose(salida99);
fclose(salida100);
    return 0;
}







//************************************************************
// FUNCTIONS
//************************************************************
void MOTION(void){

    double rp,S1,S2,S3,S4,sn,cn,sn_theta,cs_theta,Er,E_theta;
    double Br_hf,B_theta_hf,Bx_hf,By_hf,Bz_hf;
    double uxm,uym,uzm,tb,tx,ty,tz,uxr,uyr,uzr;
    double den,sx,sy,sz,uxp,uyp,uzp,fact,xx;
    float ir,jz;
    int i,j,jj1,i1;
    double psi_0=0.,factor;


   rp=sqrt(pow(x[ie],2)+pow(y[ie],2)); // radial coordinate of the particle

    if (z[ie]<0 | z[ie]>=aL | rp>r_c){
        /* no hacer calculo del movimiento*/

        if (e_out[ie]==0){
            cont_e_out++;
            e_out[ie]=1;    
        }
    }

    else{

        if (gm[ie]>=gm_max & e_gm[ie]==0){     
            cont_e_gm++;
            e_gm[ie]=1;
            }

        ir=(rp*rl)/shr;       // Index (not integer) of the particle position
        jz=(z[ie]*rl)/shz;

        i=int(ir);            // lower integer index of the particle position
        j=int(jz);

        i1=i+1;
        jj1=j+1;
        // Wheighting for magnetic field:

        S1=(jj1-jz)*(i1-ir);
        S2=(jz-j)*(i1-ir);
        S3=(jj1-jz)*(ir-i);
        S4=(jz-j)*(ir-i);
        
    	LINEAR_INTERPOLATION();

        brp=(Br[i][j]*S1+Br[i][jj1]*S2+Br[i1][j]*S3+Br[i1][jj1]*S4);
        bzp=(Bz[i][j]*S1+Bz[i][jj1]*S2+Bz[i1][j]*S3+Bz[i1][jj1]*S4);

        // Rectangular components:

        sn_theta=y[ie]/rp;
        cs_theta=x[ie]/rp;
        sn=sin(p*M_PI*z[ie]/aL);
        cn=cos(p*M_PI*z[ie]/aL);
        
        if (rp==0)
        {
         bxp=0.;
         bzp=0.;
        }
        else
        {
         bxp=brp*(x[ie]/rp);
         byp=brp*(y[ie]/rp);
        }
        // hf Magnetic field

        xx=h*rp; 
        dJ1=j0(xx)-j1(xx)/xx;

        Br_hf=2.*g0*(p*M_PI/aL)*dJ1*cn*(cs_theta*cos(tm+psi_0)+sn_theta*sin(tm+psi_0));
        B_theta_hf=-2.*g0*(p*M_PI/aL)*(j1(xx)/xx)*cn*(sn_theta*cos(tm+psi_0)-sin(tm+psi_0)*cs_theta);


        Bx_hf=Br_hf*cs_theta-B_theta_hf*sn_theta;
        By_hf=Br_hf*sn_theta+B_theta_hf*cs_theta;
        Bz_hf=2.*g0*h*j1(xx)*sn*(cs_theta*cos(tm+psi_0)+sn_theta*sin(tm+psi_0));

       // Plane Wave approximation-B
        
        bxp=bxp+Bx_hf;
        byp=byp+By_hf;
        bzp=bzp+Bz_hf;

        // hf Electric field
        // Cylindrical components

        Er=-2.*g0*(j1(xx)/xx)*sn*(cs_theta*cos(tm+psi_0)+sn_theta*sin(tm+psi_0));
        E_theta=2.*g0*dJ1*sn*(sn_theta*cos(tm+psi_0)-sin(tm+psi_0)*cs_theta);

        // Rectangular components

        Ex=Er*cs_theta-E_theta*sn_theta;
        Ey=Er*sn_theta+E_theta*cs_theta;

        // Step 1: To calculate U- 

        uxm=Ux[ie]+Ex*dt2;
        uym=Uy[ie]+Ey*dt2;
        uzm=Uz[ie];

        // Step 2: To calculate t vector

        gmn=sqrt(1.+pow(uxm,2)+pow(uym,2)+pow(uzm,2));
        tb=dt2/gmn;  

        tx=bxp*tb;
        ty=byp*tb;
        tz=bzp*tb;

        // Step 3: To calculate U'

        uxr=uxm+uym*tz-uzm*ty;
        uyr=uym+uzm*tx-uxm*tz;
        uzr=uzm+uxm*ty-uym*tx;

        // Step 4: To calculate S vector

        den=1.+pow(tx,2)+pow(ty,2)+pow(tz,2);

        sx=2.*tx/den;
        sy=2.*ty/den;
        sz=2.*tz/den;

        // Step 5:To calculate U+

        uxp=uxm+uyr*sz-uzr*sy;
        uyp=uym+uzr*sx-uxr*sz;
        uzp=uzm+uxr*sy-uyr*sx;

        // Step 6:To calculate  U n+1/2

        Ux[ie]=uxp+Ex*dt2;
        Uy[ie]=uyp+Ey*dt2;
        Uz[ie]=uzp;

        // Step 7: To calculate New coordiantes of the particle

        gm[ie]=sqrt(1+pow(Ux[ie],2)+pow(Uy[ie],2)+pow(Uz[ie],2));
        fact=dt/gm[ie];

        x[ie]=x[ie]+Ux[ie]*fact;
        y[ie]=y[ie]+Uy[ie]*fact;
        z[ie]=z[ie]+Uz[ie]*fact;

        // Calculate Phase-Shift

           PHASE_SHIFT();

        // Calculate  electromagnetic current Power transferred to the Electron Beam

        sum_energy=sum_energy+(Ex*(Ux[ie]/gm[ie])+Ey*(Uy[ie]/gm[ie]))*dt;           // Energy trasferred by TE11P mode

        rp=sqrt(pow(x[ie],2)+pow(y[ie],2));

    }

}

//**************************************************************
void PHASE_SHIFT(void)  
                       
{
    float ratio;
    int int_part;

    PHASE_VELOCITY();
    PHASE_ELECTRIC_FIELD();
    ang=ang2-ang1;     

    ratio=ang/(2.*M_PI);
    int_part=int(ratio);

    ang=(ratio-int_part)*2.*M_PI;

    if(ang<0)
    {
        ang=ang+2*M_PI;
    }
}

void PHASE_ELECTRIC_FIELD(void)
{
    if(Ex==0 & Ey>0)
    {
        ang1=M_PI/2.;
    }
    else if(Ex==0 & Ey<0)
    {
        ang1=3.*M_PI/2.;
    }
    else if(Ex>0 & Ey>=0)
    {
        ang1=atan(Ey/Ex);
    }
    else if(Ey>=0)
    {
        ang1=atan(Ey/Ex)+M_PI;
    }
    else if(Ex<0)
    {
        ang1=atan(Ey/Ex)+M_PI;
    }
    else if(Ex>0)
    {
        ang1=atan(Ey/Ex)+2.*M_PI;
    }

        ang1_new=ang1;
        ang1=ang1+M_PI;       // Because g0 is negative.
        dif1=ang1_new-ang1_old;

        if(dif1<0)
        {
            n1=n1+1;
        }

        ang1=ang1+n1*2.*M_PI;
        ang1_old=ang1_new;
}
//---------------------------
void PHASE_VELOCITY(void)
{
    if(Ux[ie]==0 & Uy[ie]>0)
    {
        ang2=M_PI/2.;
    }
    else if(Ux[ie]==0 & Uy[ie]<0)
    {
        ang2=3.*M_PI/2.;
    }
    else if(Ux[ie]>0 & Uy[ie]>=0)
    {
        ang2=atan(Uy[ie]/Ux[ie]);
    }
    else if(Uy[ie]>=0)
    {
        ang2=atan(Uy[ie]/Ux[ie])+M_PI;
    }
    else if(Ux[ie]<0)
    {
        ang2=atan(Uy[ie]/Ux[ie])+M_PI;
    }
    else if(Ux[ie]>0)
    {
        ang2=atan(Uy[ie]/Ux[ie])+2.*M_PI;
    }

        ang2_new=ang2;
        dif2=ang2_new-ang2_old;

        if(dif2<0)
        {
            n2=n2+1;
        }

        ang2=ang2+n2*2.*M_PI;
        ang2_old=ang2_new;
}


void READ_FILES(void){

    FILE *readr;
    FILE *readz;

    //printf("\n\nLeyendo Archivos de Campo.\n");

    for(int ii=1;ii<N;ii++){
      
      char nombre[32];
      char nombrefichero[64];

      strcpy(nombre, "br_component");
      sprintf(nombrefichero, "%s%d.txt",nombre,ii);
      readr=fopen(nombrefichero,"r");
      //printf("%s\n",nombrefichero );

      strcpy(nombre, "bz_component");
      sprintf(nombrefichero, "%s%d.txt",nombre,ii);
      readz=fopen(nombrefichero,"r");
      //printf("%s\n",nombrefichero );


      for(int j=0;j<npz;j++){
        for(int i=0;i<npr;i++){

          fscanf(readr,"%e",&Brr[ii][i][j]);

          fscanf(readz,"%e",&Bzr[ii][i][j]);

        }

      }

      fclose(readr);
      fclose(readz);

	
    }

    printf("\nLECTURA DE CAMPOS FINALIZADA CON ÉXITO.\n\n\n");

  }

void LINEAR_INTERPOLATION(){

    double ir;  // Índice real de la posición del ión en r.
    int i;  // índice entero de la posición del ión en r.
    double jr;  // Índice real de la posición del ión en z.
    int j;  // índice entero de la posición del ión en z.
    double ikk;	//solo son dos campos k es estatico K=1

    ri=sqrt((x[ie]*x[ie])+(y[ie]*y[ie]));
    zi=z[ie];

    ir=(ri*rl/shr);
    i=fabs(ir);

    jr=(zi*rl/shz);
    j=fabs(jr);

    pr=i;
    pz=j;

    //printf("%e \n", ik);
    double fr1[N],fr2[N],fr3[N],fr4[N];
    double fz1[N],fz2[N],fz3[N],fz4[N];



    for(int i=1;i<N;i++){

      fr1[i]=Brr[i][pr][pz];
      fr2[i]=Brr[i][pr+1][pz];
      fr3[i]=Brr[i][pr][pz+1];
      fr4[i]=Brr[i][pr+1][pz+1];

      fz1[i]=Bzr[i][pr][pz];
      fz2[i]=Bzr[i][pr+1][pz];
      fz3[i]=Bzr[i][pr][pz+1];
      fz4[i]=Bzr[i][pr+1][pz+1];

    }
    
	ikk=ik*pow(kt,-1);

    Br[pr][pz] = (fr1[2]*ikk) + (fr1[1]*(1-ikk));
    Br[pr+1][pz]= (fr2[2]*ikk) + (fr2[1]*(1-ikk));
    Br[pr][pz+1]= (fr3[2]*ikk) + (fr3[1]*(1-ikk));
    Br[pr+1][pz+1]= (fr4[2]*ikk) + (fr4[1]*(1-ikk));

    Bz[pr][pz]= (fz1[2]*ikk) + (fz1[1]*(1-ikk));
    Bz[pr+1][pz]= (fz2[2]*ikk) + (fz2[1]*(1-ikk));
    Bz[pr][pz+1]= (fz3[2]*ikk) + (fz3[1]*(1-ikk));
    Bz[pr+1][pz+1]= (fz4[2]*ikk) + (fz4[1]*(1-ikk));

  }

void INJECT_ELECTRON(void)              // (1) Random Distribution
{
 //  Inyectar un haz de electrones
    double rdm1,rdm2,rdm3,rdm4;
    double xp,yp,zp,rp,Lz0=vz0*dt;//Lz0=ki*vz0*dt;
    double R_ball=0.5*r_c;       
    int seed = time (NULL); srand (seed);
    double vx,vy,vz, gama;
    rp=0.;
    check=0;
   
      for(int i=0;i<Ne;i++) 
        {      

            while(1){
                rdm1=rand();
                rdm1=rdm1/RAND_MAX;
                rdm2=rand();
                rdm2=rdm2/RAND_MAX;
                rdm3=rand();
                rdm3=rdm3/RAND_MAX;

                xp=(2.*rdm1-1.0)*R_ball;
                yp=(2.*rdm2-1.0)*R_ball;
                zp=(2.*rdm3-1.0)*R_ball;

                rp=sqrt(pow(xp,2)+pow(yp,2)+pow(zp,2));
                zp=aL/2 + zp;

                if(rp<R_ball)
                {
                    check=check+1;


                    // Electrons
                        x[i]=xp;
                        y[i]=yp;
                        z[i]=zp;

                        x[99]=x0;
                        y[99]=y00;
                        z[99]=z0;

                        Ux[i]=distrib_vel()*(vte/c); //la velocidad que se coloca en z ahora esta en x 
                        Uy[i]=distrib_vel()*(vte/c);
                        Uz[i]=distrib_vel()*(vte/c);
                        Ux[99]=0.0001;
                        Uy[99]=0.0001;
                        Uz[99]=0.;
                    
                    break;
                }

            }
        }     

}

double distrib_vel (void) // función para generar distribución maxwelliana de la velocidad de las particulas 
                              // (Ver pág. 291 Computational Physics Fitzpatrick: Distribution functions--> Rejection Method)
{

 double sigma=1;                  // sigma=vt/vt (normalyzed)
 double vmin=-4*sigma;           // Rapidez mínima  
 double vmax= 4*sigma;           // Rapidez máxima
 double n0_1D=pow(Ne,1/3.)*100;   // Concentración 1D (1/m)
 double ni_1D=n0_1D*lambda_Debye; // Concentración 1D (normalizada)
 double fmax=ni_1D/sqrt(2*M_PI);  // Densidad de propabilidad máxima
 double v,vx,f,f_random;


//  Inizializar generador de números aleatorios 
 
static int flag = 0;

if (flag == 0)
  {
    int seed = time (NULL);
    srand (seed);
    flag = 1;
  }

  v=vmin+(vmax-vmin)*double(rand())/RAND_MAX;                        // Calcular valor aleatorio de v uniformemente distribuido en el rango [vmin,vmax]
  f =ni_1D*(exp(-pow(v,2)/2.))/sqrt(2*M_PI);                        // Evaluar f en el v calculado.                     
  f_random = fmax*double(rand())/RAND_MAX;                           // Calcular valor aleatorio de f uniformemente distribuido en el rango [0,fmax]
  if (f_random > f) return distrib_vel();                          // If f_random > f rechazar el valor y recalcular
  else return v;
}


void PRINT_FIELD(void){

 for(int j=0;j<npz;j++){
    fprintf(out11,"%d	%e	%e\n",j,Bzr[1][0][j],Bzr[2][0][j]); // imprime campo magnetico en una curva de nivel 
       for(int ix=0;ix<=2*npr-2;ix++){
		

            if (ix<=npr-1)
            {
                Bxz_2=pow(Brr[1][npr-1-ix][j],2)+pow(Bzr[1][npr-1-ix][j],2);
                fprintf(out2,"%e      ",sqrt(Bxz_2)); // (kV)

                Bxz_2=pow(Brr[2][npr-1-ix][j],2)+pow(Bzr[2][npr-1-ix][j],2);
                fprintf(out3,"%e      ",sqrt(Bxz_2)); // (kV)
            }
            else
            {
                Bxz_2=pow(Brr[1][ix-npr+1][j],2)+pow(Bzr[1][ix-npr+1][j],2);
                fprintf(out2,"%e      ",sqrt(Bxz_2));

                Bxz_2=pow(Brr[2][ix-npr+1][j],2)+pow(Bzr[2][ix-npr+1][j],2);
                fprintf(out3,"%e      ",sqrt(Bxz_2));   
            }
        }
        fprintf(out2,"\n");
        fprintf(out3,"\n");
    }
    fclose(out2);
    fclose(out3);
    fclose(out11);	
}

void PRINT_POSITION(void){


/*
   for (int i = 0; i < Ne; i++){
       sprintf(buffer, "%02d.txt", i+1); // Ahora tenemos en buffer = "holaXX.txt"
       fopen(Posiciones/,buffer, "w");
    }*/
fprintf(salida1,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[0]*rl*100.,y[0]*rl*100.,z[0]*rl*100.,gm[0]);
fprintf(salida2,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[1]*rl*100.,y[1]*rl*100.,z[1]*rl*100.,gm[1]);
fprintf(salida3,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[2]*rl*100.,y[2]*rl*100.,z[2]*rl*100.,gm[2]);
fprintf(salida4,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[3]*rl*100.,y[3]*rl*100.,z[3]*rl*100.,gm[3]);
fprintf(salida5,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[4]*rl*100.,y[4]*rl*100.,z[4]*rl*100.,gm[4]);
fprintf(salida6,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[5]*rl*100.,y[5]*rl*100.,z[5]*rl*100.,gm[5]);
fprintf(salida7,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[6]*rl*100.,y[6]*rl*100.,z[6]*rl*100.,gm[6]);
fprintf(salida8,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[7]*rl*100.,y[7]*rl*100.,z[7]*rl*100.,gm[7]);
fprintf(salida9,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[8]*rl*100.,y[8]*rl*100.,z[8]*rl*100.,gm[8]);
fprintf(salida10,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[9]*rl*100.,y[9]*rl*100.,z[9]*rl*100.,gm[9]);
fprintf(salida11,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[10]*rl*100.,y[10]*rl*100.,z[10]*rl*100.,gm[10]);
fprintf(salida12,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[11]*rl*100.,y[11]*rl*100.,z[11]*rl*100.,gm[11]);
fprintf(salida13,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[12]*rl*100.,y[12]*rl*100.,z[12]*rl*100.,gm[12]);
fprintf(salida14,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[13]*rl*100.,y[13]*rl*100.,z[13]*rl*100.,gm[13]);
fprintf(salida15,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[14]*rl*100.,y[14]*rl*100.,z[14]*rl*100.,gm[14]);
fprintf(salida16,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[15]*rl*100.,y[15]*rl*100.,z[15]*rl*100.,gm[15]);
fprintf(salida17,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[16]*rl*100.,y[16]*rl*100.,z[16]*rl*100.,gm[16]);
fprintf(salida18,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[17]*rl*100.,y[17]*rl*100.,z[17]*rl*100.,gm[17]);
fprintf(salida19,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[18]*rl*100.,y[18]*rl*100.,z[18]*rl*100.,gm[18]);
fprintf(salida20,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[19]*rl*100.,y[19]*rl*100.,z[19]*rl*100.,gm[19]);
fprintf(salida21,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[20]*rl*100.,y[20]*rl*100.,z[20]*rl*100.,gm[20]);
fprintf(salida22,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[21]*rl*100.,y[21]*rl*100.,z[21]*rl*100.,gm[21]);
fprintf(salida23,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[22]*rl*100.,y[22]*rl*100.,z[22]*rl*100.,gm[22]);
fprintf(salida24,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[23]*rl*100.,y[23]*rl*100.,z[23]*rl*100.,gm[23]);
fprintf(salida25,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[24]*rl*100.,y[24]*rl*100.,z[24]*rl*100.,gm[24]);
fprintf(salida26,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[25]*rl*100.,y[25]*rl*100.,z[25]*rl*100.,gm[25]);
fprintf(salida27,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[26]*rl*100.,y[26]*rl*100.,z[26]*rl*100.,gm[26]);
fprintf(salida28,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[27]*rl*100.,y[27]*rl*100.,z[27]*rl*100.,gm[27]);
fprintf(salida29,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[28]*rl*100.,y[28]*rl*100.,z[28]*rl*100.,gm[28]);
fprintf(salida30,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[29]*rl*100.,y[29]*rl*100.,z[29]*rl*100.,gm[29]);
fprintf(salida31,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[30]*rl*100.,y[30]*rl*100.,z[30]*rl*100.,gm[30]);
fprintf(salida32,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[31]*rl*100.,y[31]*rl*100.,z[31]*rl*100.,gm[31]);
fprintf(salida33,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[32]*rl*100.,y[32]*rl*100.,z[32]*rl*100.,gm[32]);
fprintf(salida34,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[33]*rl*100.,y[33]*rl*100.,z[33]*rl*100.,gm[33]);
fprintf(salida35,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[34]*rl*100.,y[34]*rl*100.,z[34]*rl*100.,gm[34]);
fprintf(salida36,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[35]*rl*100.,y[35]*rl*100.,z[35]*rl*100.,gm[35]);
fprintf(salida37,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[36]*rl*100.,y[36]*rl*100.,z[36]*rl*100.,gm[36]);
fprintf(salida38,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[37]*rl*100.,y[37]*rl*100.,z[37]*rl*100.,gm[37]);
fprintf(salida39,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[38]*rl*100.,y[38]*rl*100.,z[38]*rl*100.,gm[38]);
fprintf(salida40,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[39]*rl*100.,y[39]*rl*100.,z[39]*rl*100.,gm[39]);
fprintf(salida41,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[40]*rl*100.,y[40]*rl*100.,z[40]*rl*100.,gm[40]);
fprintf(salida42,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[41]*rl*100.,y[41]*rl*100.,z[41]*rl*100.,gm[41]);
fprintf(salida43,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[42]*rl*100.,y[42]*rl*100.,z[42]*rl*100.,gm[42]);
fprintf(salida44,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[43]*rl*100.,y[43]*rl*100.,z[43]*rl*100.,gm[43]);
fprintf(salida45,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[44]*rl*100.,y[44]*rl*100.,z[44]*rl*100.,gm[44]);
fprintf(salida46,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[45]*rl*100.,y[45]*rl*100.,z[45]*rl*100.,gm[45]);
fprintf(salida47,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[46]*rl*100.,y[46]*rl*100.,z[46]*rl*100.,gm[46]);
fprintf(salida48,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[47]*rl*100.,y[47]*rl*100.,z[47]*rl*100.,gm[47]);
fprintf(salida49,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[48]*rl*100.,y[48]*rl*100.,z[48]*rl*100.,gm[48]);
fprintf(salida50,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[49]*rl*100.,y[49]*rl*100.,z[49]*rl*100.,gm[49]);
fprintf(salida51,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[50]*rl*100.,y[50]*rl*100.,z[50]*rl*100.,gm[50]);
fprintf(salida52,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[51]*rl*100.,y[51]*rl*100.,z[51]*rl*100.,gm[51]);
fprintf(salida53,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[52]*rl*100.,y[52]*rl*100.,z[52]*rl*100.,gm[52]);
fprintf(salida54,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[53]*rl*100.,y[53]*rl*100.,z[53]*rl*100.,gm[53]);
fprintf(salida55,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[54]*rl*100.,y[54]*rl*100.,z[54]*rl*100.,gm[54]);
fprintf(salida56,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[55]*rl*100.,y[55]*rl*100.,z[55]*rl*100.,gm[55]);
fprintf(salida57,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[56]*rl*100.,y[56]*rl*100.,z[56]*rl*100.,gm[56]);
fprintf(salida58,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[57]*rl*100.,y[57]*rl*100.,z[57]*rl*100.,gm[57]);
fprintf(salida59,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[58]*rl*100.,y[58]*rl*100.,z[58]*rl*100.,gm[58]);
fprintf(salida60,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[59]*rl*100.,y[59]*rl*100.,z[59]*rl*100.,gm[59]);
fprintf(salida61,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[60]*rl*100.,y[60]*rl*100.,z[60]*rl*100.,gm[60]);
fprintf(salida62,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[61]*rl*100.,y[61]*rl*100.,z[61]*rl*100.,gm[61]);
fprintf(salida63,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[62]*rl*100.,y[62]*rl*100.,z[62]*rl*100.,gm[62]);
fprintf(salida64,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[63]*rl*100.,y[63]*rl*100.,z[63]*rl*100.,gm[63]);
fprintf(salida65,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[64]*rl*100.,y[64]*rl*100.,z[64]*rl*100.,gm[64]);
fprintf(salida66,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[65]*rl*100.,y[65]*rl*100.,z[65]*rl*100.,gm[65]);
fprintf(salida67,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[66]*rl*100.,y[66]*rl*100.,z[66]*rl*100.,gm[66]);
fprintf(salida68,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[67]*rl*100.,y[67]*rl*100.,z[67]*rl*100.,gm[67]);
fprintf(salida69,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[68]*rl*100.,y[68]*rl*100.,z[68]*rl*100.,gm[68]);
fprintf(salida70,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[69]*rl*100.,y[69]*rl*100.,z[69]*rl*100.,gm[69]);
fprintf(salida71,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[70]*rl*100.,y[70]*rl*100.,z[70]*rl*100.,gm[70]);
fprintf(salida72,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[71]*rl*100.,y[71]*rl*100.,z[71]*rl*100.,gm[71]);
fprintf(salida73,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[72]*rl*100.,y[72]*rl*100.,z[72]*rl*100.,gm[72]);
fprintf(salida74,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[73]*rl*100.,y[73]*rl*100.,z[73]*rl*100.,gm[73]);
fprintf(salida75,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[74]*rl*100.,y[74]*rl*100.,z[74]*rl*100.,gm[74]);
fprintf(salida76,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[75]*rl*100.,y[75]*rl*100.,z[75]*rl*100.,gm[75]);
fprintf(salida77,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[76]*rl*100.,y[76]*rl*100.,z[76]*rl*100.,gm[76]);
fprintf(salida78,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[77]*rl*100.,y[77]*rl*100.,z[77]*rl*100.,gm[77]);
fprintf(salida79,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[78]*rl*100.,y[78]*rl*100.,z[78]*rl*100.,gm[78]);
fprintf(salida80,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[79]*rl*100.,y[79]*rl*100.,z[79]*rl*100.,gm[79]);
fprintf(salida81,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[80]*rl*100.,y[80]*rl*100.,z[80]*rl*100.,gm[80]);
fprintf(salida82,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[81]*rl*100.,y[81]*rl*100.,z[81]*rl*100.,gm[81]);
fprintf(salida83,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[82]*rl*100.,y[82]*rl*100.,z[82]*rl*100.,gm[82]);
fprintf(salida84,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[83]*rl*100.,y[83]*rl*100.,z[83]*rl*100.,gm[83]);
fprintf(salida85,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[84]*rl*100.,y[84]*rl*100.,z[84]*rl*100.,gm[84]);
fprintf(salida86,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[85]*rl*100.,y[85]*rl*100.,z[85]*rl*100.,gm[85]);
fprintf(salida87,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[86]*rl*100.,y[86]*rl*100.,z[86]*rl*100.,gm[86]);
fprintf(salida88,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[87]*rl*100.,y[87]*rl*100.,z[87]*rl*100.,gm[87]);
fprintf(salida89,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[88]*rl*100.,y[88]*rl*100.,z[88]*rl*100.,gm[88]);
fprintf(salida90,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[89]*rl*100.,y[89]*rl*100.,z[89]*rl*100.,gm[89]);
fprintf(salida91,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[90]*rl*100.,y[90]*rl*100.,z[90]*rl*100.,gm[90]);
fprintf(salida92,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[91]*rl*100.,y[91]*rl*100.,z[91]*rl*100.,gm[91]);
fprintf(salida93,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[92]*rl*100.,y[92]*rl*100.,z[92]*rl*100.,gm[92]);
fprintf(salida94,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[93]*rl*100.,y[93]*rl*100.,z[93]*rl*100.,gm[93]);
fprintf(salida95,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[94]*rl*100.,y[94]*rl*100.,z[94]*rl*100.,gm[94]);
fprintf(salida96,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[95]*rl*100.,y[95]*rl*100.,z[95]*rl*100.,gm[95]);
fprintf(salida97,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[96]*rl*100.,y[96]*rl*100.,z[96]*rl*100.,gm[96]);
fprintf(salida98,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[97]*rl*100.,y[97]*rl*100.,z[97]*rl*100.,gm[97]);
fprintf(salida99,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[98]*rl*100.,y[98]*rl*100.,z[98]*rl*100.,gm[98]);
fprintf(salida100,"%e %e %e %e %e\n",ik*pow(kt,-1)*4,56,x[99]*rl*100.,y[99]*rl*100.,z[99]*rl*100.,gm[99]);

    if (ik==0){
        for (int i = 0; i < Ne; ++i){
            if (e_out[i]==1){} //no hacer nada
            else{    
             fprintf(out4,"%e ",x[i]*rl*100.);
             fprintf(out4,"%e ",y[i]*rl*100.);
             fprintf(out4,"%e ",z[i]*rl*100.);
	     fprintf(out4,"%e ",gm[i]);
             fprintf(out4,"%e\n",(gm[i]-1.)*511.);
            }
        }
        fclose(out4);  
    } 
    
    if (ik==kt/4){
        for (int i = 0; i < Ne; ++i){
            if (e_out[i]==1){} //no hacer nada
            else{ 
             fprintf(out5,"%e ",x[i]*rl*100.);
             fprintf(out5,"%e ",y[i]*rl*100.);
             fprintf(out5,"%e ",z[i]*rl*100.);
	     fprintf(out5,"%e ",gm[i]);
             fprintf(out5,"%e\n",(gm[i]-1.)*511.);
            }
        }
        fclose(out5);
    }     

    if (ik==kt/2){
        for (int i = 0; i < Ne; ++i){
            if (e_out[i]==1){} //no hacer nada
            else{             
             fprintf(out6,"%e ",x[i]*rl*100.);
             fprintf(out6,"%e ",y[i]*rl*100.);
             fprintf(out6,"%e ",z[i]*rl*100.);
	     fprintf(out6,"%e ",gm[i]);
             fprintf(out6,"%e\n",(gm[i]-1.)*511.);
            }             
        }
        fclose(out6);
    }  

    if (ik==kt/4*3){
        for (int i = 0; i < Ne; ++i){
            if (e_out[i]==1){} //no hacer nada
            else{             
             fprintf(out7,"%e ",x[i]*rl*100.);
             fprintf(out7,"%e ",y[i]*rl*100.);
             fprintf(out7,"%e ",z[i]*rl*100.);
	     fprintf(out7,"%e ",gm[i]);
             fprintf(out7,"%e\n",(gm[i]-1.)*511.);
            }             
        }
        fclose(out7);
    }  

    if (ik==kt-1){

	 gmax=1.0;

        for (int i = 0; i < Ne; ++i){
            if (e_out[i]==1){} //no hacer nada
            else{             
					
             fprintf(out8,"%e ",x[i]*rl*100.);
             fprintf(out8,"%e ",y[i]*rl*100.);
             fprintf(out8,"%e ",z[i]*rl*100.);
	     	 fprintf(out8,"%e ",gm[i]);
             fprintf(out8,"%e\n",(gm[i]-1.)*511.);
            }             
        }
        fclose(out8);
    } 
}
