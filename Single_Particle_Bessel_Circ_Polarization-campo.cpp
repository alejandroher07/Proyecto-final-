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


#define Ne 100      //numero de electrones           
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
double m_gm=16;	//maximo valor de gamma alcanzado en el sistema (para realizar el espectro energetico)

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




FILE *out,*out2,*out3,*out4,*out5,*out6,*out7,*out8,*out9,*out10,*out11;
int main()
{
    
    out=fopen("Results.txt","w");
    out2=fopen("B_rz_i.txt","w");
    out3=fopen("B_rz_f.txt","w");
    out4=fopen("Posiciones_0.txt","w");
    out5=fopen("Posiciones_25.txt","w");
    out6=fopen("Posiciones_50.txt","w");
    out7=fopen("Posiciones_75.txt","w");
    out8=fopen("Posiciones_100.txt","w");
    out9=fopen("Espectro_energetico_inicial.txt","w");
    out10=fopen("Espectro_energetico_final.txt","w");
    out11=fopen("B_z_i_f.txt","w");

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
		gmax=gmax+(m_gm/npe);

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


		if(cont_e_gm>=10){
			ref=ik;
			for(int i=ref;i<2*ref;i++){
				ik=ik-1;
				tm=ik*dt;          // current time
				kk=ik-kp;          // Index for electric field
				for (ie = 0; ie < Ne; ++ie){

					if(e_gm[ie]==1){
						printf("x=%e y=%e z=%e	ik=%d\n",x[ie]*rl*100,y[ie]*rl*100,z[ie]*rl*100,ik);
						MOTION();
						}

					}
 				   if(ik==kp+m-1) kp=kp+m;
				}
		
		exit(1);
		}

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

                        Ux[i]=distrib_vel()*(vte/c); //la velocidad que se coloca en z ahora esta en x 
                        Uy[i]=distrib_vel()*(vte/c);
                        Uz[i]=distrib_vel()*(vte/c);
                    
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















    if (ik==0){
        for (int i = 0; i < Ne; ++i){
            if (e_out[i]==1){} //no hacer nada
            else{    
	     for(int ii=1;ii<npe;ii++)if(gm[i]>g_max[ii-1] & gm[i]<g_max[ii])ene_cont[ii]=ene_cont[ii]+1; 	//ciclo para hacer espectro energetico	
		
             fprintf(out4,"%e    ",x[i]*rl*100.);
             fprintf(out4,"%e    ",y[i]*rl*100.);
             fprintf(out4,"%e    ",z[i]*rl*100.);
	     fprintf(out4,"%e    ",gm[i]);
             fprintf(out4,"%e    \n",(gm[i]-1.)*511.);
            }
        }
        fclose(out4);  
	for(int ii=1;ii<npe;++ii)fprintf(out9,"%e    %e    %e    \n",g_max[ii],(g_max[ii]-1.)*511,ene_cont[ii]/npe);
	fclose(out9);
    } 
    
    if (ik==kt/4){
        for (int i = 0; i < Ne; ++i){
            if (e_out[i]==1){} //no hacer nada
            else{ 
             fprintf(out5,"%e    ",x[i]*rl*100.);
             fprintf(out5,"%e    ",y[i]*rl*100.);
             fprintf(out5,"%e    ",z[i]*rl*100.);
	     fprintf(out5,"%e    ",gm[i]);
             fprintf(out5,"%e    \n",(gm[i]-1.)*511.);
            }
        }
        fclose(out5);
    }     

    if (ik==kt/2){
        for (int i = 0; i < Ne; ++i){
            if (e_out[i]==1){} //no hacer nada
            else{             
             fprintf(out6,"%e    ",x[i]*rl*100.);
             fprintf(out6,"%e    ",y[i]*rl*100.);
             fprintf(out6,"%e    ",z[i]*rl*100.);
	     fprintf(out6,"%e    ",gm[i]);
             fprintf(out6,"%e    \n",(gm[i]-1.)*511.);
            }             
        }
        fclose(out6);
    }  

    if (ik==kt/4*3){
        for (int i = 0; i < Ne; ++i){
            if (e_out[i]==1){} //no hacer nada
            else{             
             fprintf(out7,"%e    ",x[i]*rl*100.);
             fprintf(out7,"%e    ",y[i]*rl*100.);
             fprintf(out7,"%e    ",z[i]*rl*100.);
	     fprintf(out7,"%e    ",gm[i]);
             fprintf(out7,"%e    \n",(gm[i]-1.)*511.);
            }             
        }
        fclose(out7);
    }  

    if (ik==kt-1){

	 gmax=1.0;
         for (int i = 0; i < npe; ++i){  
		ene_cont[i]=0; //inicializacion para conteo energetico
		g_max[i]=gmax;
		gmax=gmax+(m_gm/npe);
         }

        for (int i = 0; i < Ne; ++i){
            if (e_out[i]==1){} //no hacer nada
            else{             

	     for(int ii=1;ii<npe;ii++)if(gm[i]>g_max[ii-1] & gm[i]<g_max[ii])ene_cont[ii]=ene_cont[ii]+1; 	//ciclo para hacer espectro energetico		
					
             fprintf(out8,"%e    ",x[i]*rl*100.);
             fprintf(out8,"%e    ",y[i]*rl*100.);
             fprintf(out8,"%e    ",z[i]*rl*100.);
	     fprintf(out8,"%e    ",gm[i]);
             fprintf(out8,"%e    \n",(gm[i]-1.)*511.);
            }             
        }
        fclose(out8);
	for(int ii=1;ii<npe;++ii)fprintf(out10,"%e    %e    %e    \n",g_max[ii],(g_max[ii]-1.)*511,ene_cont[ii]/npe);
	fclose(out10);
    } 
}
