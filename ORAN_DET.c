#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "/Users/javo/libreria/aleatorios.h"
#include "/Users/javo/libreria/SpecialFunctions.h"
#include "/Users/javo/libreria/ranlib.h"
#include "/Users/javo/libreria/rnglib.h"
//#include "/Users/javo/libreria/libmul.h"
#include "/Users/javo/libreria/ode.c"


#include "def_oran.c"

#define m 13 // numero de ecuaciones

void sistema(double t, double v[], double dv[]); 
/* sourv exp */
double suv_exp(double tt,double s, double rate);
/*Fbar - supervivencia*/
double Fbar(double tt,double s, double k, double theta);

/* variables globales */

double 
	beta_day_theta_0, fR, mu_Dry,
	m_E_C_G, mu_Wet,
	m_L, mu_L, C_L,
	m_P, mu_P,
	m_M, mu_M,
	mu_V,
	GAMA, b_theta_pV,
	EV, sigma_V,
	b_theta_pH,
	sigma_H,
	gama, deltaI;

	


/* Comienza el programa */
int main(int argc,char *argv[])
{
	FILE *archivo;
    FILE *vector_host;
    FILE *casos;
    FILE *mosquitos;
    
    float pop_Vect  = 0.;
    mosquitos		=   fopen("mosquitos_det.txt","w");
    
    vector_host 	= 	fopen("v_h_det.txt","w");
    casos 			= 	fopen("casos_det.txt","w");
    
    float v_h_week = 0.;
    
    /*year, moth, week, day*/
    int 
    day		= 0,
    week 	= 0;/*,
    year	= BEGIN_YEAR,
    month	= 365;*/
	
	/*count */
    int i,j;
    
    double t=0., dt =  Dt; //one day
    
    /* variables ambientales */
    float Tmin[DAYS],
          Tmax[DAYS],
          Tmean[DAYS],
          Tmedian[DAYS],
          Rain[DAYS],
          Hum_mean[DAYS],
          Hum_median[DAYS];
          
    float Tm=0.,RH=0.,R=0.,Tdeath;
    
    double T_1, media_VE_1, var_VE_1, theta_VE_1,
		k_VE_1, mu_V_1, sigma_V_1, count_s,
		T_2, media_VE_2, sigma_V_2, mu_V_2 ;
          
    double v[m],dv[m];
          
    double integral_1=0.,integral_2=0.,
		integral_3=0., integral_4=0.;
    
    double G_T[DAYS],mu_U_T[DAYS],sigma_U_T[DAYS],U_T[DAYS];
    
    double pob, media_VE,var_VE,theta_VE,k_VE;
    
    float case_day, case_week;
    
    float hogares 	= 0.,//17650.,
			poblacion	= 0.;//76070.;
			
	float Host[RAD_CENSAL],Hog[RAD_CENSAL];
    
    
    /*eviromental values*/
      archivo = fopen("Oran_2007_2017_medio.txt","r");
      for(i=0;i<DAYS;i++) {fscanf(archivo,"%f\t%f\t%f\t%f\t%f\t%f\t%f\n",&Tmin[i],&Tmax[i],&Tmean[i],&Tmedian[i],&Rain[i],&Hum_mean[i],&Hum_median[i]);}
	  fclose(archivo);
	  
	 /*Hogares por radio censal*/
	  archivo = fopen("hogares_oran.txt","r");
      for(  j = 0; j < RAD_CENSAL; j++  ){
        fscanf( archivo, "%f\n", &Hog[j]);
        hogares = hogares + Hog[j];
        }
	  fclose(archivo);
	  
	  /*array census radius used ; as array*/
	  archivo = fopen("censo_oran.txt","r");
      for(  j = 0; j < RAD_CENSAL; j++  ){
        fscanf( archivo, "%f\n", &Host[j] );
        poblacion = poblacion + Host[j];
      }
      fclose(archivo);
	  
	  
	  /* PARA EL MODELO DE CASOS IMPORTADOS */
	  
	  float rate_week_imp[WEEKS], imported[DAYS],importados[DAYS];
	  
	  /*infection rate import per week 0~1.*/
      archivo = fopen("serie_bolivia.txt","r");
      for(i=0;i<WEEKS;i++){fscanf(archivo,"%f\n",&rate_week_imp[i]);}
      fclose(archivo);
      /*imported people*/
	  archivo = fopen("aguas.txt","r");
      for(i=0;i<DAYS;i++) {fscanf(archivo,"%f\n",&imported[i]);}
	  fclose(archivo);
	 int se = 0;//semana epi
	 //int semana = 7;
	 float inf_import,aux0=0.,inf_import_bol;
     archivo = fopen("CI_ORAN_DET.txt","w");
     
     static long int seed0;	
	 seed0 = -atol(argv[1]); //usar siempre 12

	  for( i = 0; i < DAYS; i++ ){
          day = day + 1;
		  inf_import = (float)(RATE_IMPORT*(rate_week_imp[se]*imported[i])/POPULATION);
  //       if (i == DAYS - 1) { printf("%d\t%d\t%.0f\n",corr,i,importados[i]);
		  inf_import_bol = poidev( dt*inf_import, &seed0 );
          aux0 = aux0 + inf_import_bol;
          importados[i] = inf_import_bol;
          if ( day == 7 ){
			  se++;
              fprintf(archivo,"%d\t%.2f\n",se,aux0);
              aux0 = 0.;
			  day = 0;
			  }
          //if (inf_import>0){printf("%d\t%.2f\n",se,inf_import);}
        }
    fclose(archivo);
	  
	  /* lectura cond iniciales */
    float ED0[RAD_CENSAL],EW0[RAD_CENSAL],Larv0[RAD_CENSAL],Pupa0[RAD_CENSAL],MOSCO0[RAD_CENSAL],H_t0[RAD_CENSAL];
    /* initial conditions censal radius*/
	  archivo = fopen("bichos.txt","r");
      for(j=0;j<RAD_CENSAL;j++){fscanf(archivo,"%f\t%f\t%f\t%f\t%f\t%f\n",&ED0[j],&EW0[j],&Larv0[j],&Pupa0[j],&MOSCO0[j],&H_t0[j]);} //radios censales subidos
      fclose(archivo);
	  
			
			/* variables auxiliares*/
			
	  float KL, Larv, H_t = 24.;

    
		/*indice*/
     for ( i = 0; i < DAYS; i++){
		 
		 G_T[i] 		= 	0.;
		 U_T[i] 		= 	0.;
		 mu_U_T[i] 		= 	0.;
		 sigma_U_T[i] 	= 	0.;
		 
	 }
    
	
	int count = 0;
	
	/* Condiciones iniciales */
	
	t 		=	1.;
	
	day		=	0.;
	
	week	=	1;
	
	H_t = 24.;
	
	v[0] = 0.; //E_D 		
	v[1] = 0.; //E_W		
	v[2] = 0.; //L		
	v[3] = 0.; //P		
	v[4] = 0.; //M		
	v[5] = 0.; //V		
	v[6] = 0.; //V_S		
	v[7] = 0.; //V_E		
	v[8] = 0.; //V_I		
	v[9] = ALPHA*poblacion; //H_S		
	v[10] = 0.; //H_E		
	v[11] = 0.; //H_I		
	v[12] = poblacion - v[9] - v[10] - v[11]; //H_R		

	for(j=0;j<RAD_CENSAL;j++){
		
		v[0] 	= 	v[0] + ED0[j]; //E_D 		
		v[1] 	= 	v[0] + EW0[j]; //E_W		
		v[2] 	= 	v[0] + Larv0[j]; //L		
		v[3] 	= 	v[0] + Pupa0[j]; //P		
		v[4] 	= 	v[0] + MOSCO0[j]; //M	
	}
	
	
	

	
	while ( t < DAYS ){
		
		/* calculo EV */
		
		Tm			=	Tmean[count];
		
		media_VE	=	1+(0.1216*Tm*Tm - 8.66*Tm + 154.79);
		
		var_VE		=	1+(0.1728*Tm*Tm - 12.36*Tm + 230.62);
		
		sigma_V		=	1./media_VE;
		
		k_VE		=	( media_VE*media_VE)/var_VE;
		
		theta_VE	=	var_VE / media_VE;
		
		mu_V		=	muerte_V(Tm)*MU_MOSQUITA_ADULTA;
		
		pob			=	v[9] + v[10] + v[11] +v[12];
		
		G_T[count]	=  b_theta_pV * v[6] * (v[11]/pob);
		
		if ( count > 3){
			
			for ( i = 1; i < count; i++){
				
				T_1				=	Tmean[i];
				
				media_VE_1		=	0.1216*T_1*T_1 - 8.66*T_1 + 154.79;
		
				var_VE_1		=	0.1728*T_1*T_1 - 12.36*T_1 + 230.62;
		
				sigma_V_1		=	1./media_VE_1;
		
				k_VE_1			=	( media_VE_1*media_VE_1)/var_VE_1;
		
				theta_VE_1		=	var_VE_1 / media_VE_1;
		
				mu_V_1			=	muerte_V((float)T_1)*MU_MOSQUITA_ADULTA;
				
				count_s			= 	(double)count;
				
				sigma_U_T[i] 	= 	sigma_V_1*Fbar(t, count_s , k_VE_1, theta_VE_1)*suv_exp(t, count_s, mu_V_1);
				
				mu_U_T[i] 		= 	mu_V_1*Fbar(t, count_s , k_VE, theta_VE)*suv_exp(t, count_s, mu_V_1);
				
				U_T[i] 			= 	Fbar(t, count_s , k_VE, theta_VE)*suv_exp(t, count_s, mu_V_1);
				
			}
			
			integral_1	=	0.;
			integral_2	=	0.;
			integral_3	=	0.;
			integral_4	=	0.;
			
			for ( i = 1; i < count; i++){
				
				T_1				=	Tmean[i-1];
				
				media_VE_1		=	0.1216*T_1*T_1 - 8.66*T_1 + 154.79;
				
				sigma_V_1		=	1./media_VE_1;
				
				mu_V_1			=	muerte_V((float)T_1)*MU_MOSQUITA_ADULTA;
				
				T_2				=	Tmean[i];
				
				media_VE_2		=	0.1216*T_2*T_2 - 8.66*T_2 + 154.79;
				
				sigma_V_2		=	1./media_VE_2;
				
				mu_V_2			=	muerte_V((float)T_2)*MU_MOSQUITA_ADULTA;
				
			
				integral_1		=	integral_1 + 0.5*( G_T[i-1]*U_T[i-1] + G_T[i]*U_T[i] );
				if ( integral_1 < 0.){ integral_1 = 0.;}
			
				integral_2		=	integral_2 + 0.5*( sigma_V_1*G_T[i-1]*U_T[i-1] + sigma_V_2*G_T[i]*U_T[i] );
				if ( integral_2 < 0.){ integral_2 = 0.;}
			
				integral_3		=	integral_3 + 0.5*( G_T[i-1]*U_T[i-1] + G_T[i]*U_T[i] );
				if ( integral_3 < 0.){ integral_3 = 0.;}
				
				integral_4		=	integral_4 + 0.5*( mu_V_1*G_T[i-1]*U_T[i-1] + mu_V_1*G_T[i]*U_T[i] );
				if ( integral_4 < 0.){ integral_4 = 0.;}
				//printf ("%d\t%.2f\t%.2f\t%.2f\t%.2f\n",count,integral_1,integral_2,integral_3,integral_4);
			
			}
		}
		else {
			
			integral_1	=	0.;
			
			integral_2	=	0.;
			
			integral_3	=	0.;
			
			integral_4	=	0.;
			
		}
		
		/* valores metereologicos*/
		
		Tm			=	Tmean[count];
		R			=	Rain[count];
		RH			=	Hum_mean[count];
		Tdeath		=	Tmin[count];
		
		
		
		/* modelo poblacion aedes */
		
		beta_day_theta_0	=	beta_day*theta_T(Tm);
		
		fR					=	egg_wet(R);
		
		mu_Dry				=	1./EGG_LIFE;
		
		Larv		=	(float)v[2];
		
		H_t 		= 	H_t_1(R, RH, H_t, Tm ); 
		
		KL			=	hogares*( Kmax*H_t/Hmax ) + 1.0;
		
		m_E_C_G		=	0.24*rate_mx(Tm, 10798.,100000.,14184.)*C_Gillet(Larv, KL) ;
		
		mu_Wet		=	1./EGG_LIFE_wet;
		
		m_L		=	0.2088*rate_mx(Tm, 26018.,55990.,304.6);
		
		if( Tdeath < 13.4 ){  m_L = 0.; } //larvario
		
		mu_L	=	0.01 + 0.9725*exp(- (Tm - 4.85)/2.7035);
		
		C_L		=	1.5*(Larv/KL);
		
		m_P			=	0.384*rate_mx(Tm, 14931.,-472379.,148.);
		
		mu_P		=	0.01 + 0.9725*exp(- (Tm - 4.85)/2.7035);
		
		m_M		=	MADURACION_MOSQUITO;
		
		mu_M	=	MU_MOSQUITO_JOVEN;
		
		mu_V	=	muerte_V(Tm)*MU_MOSQUITA_ADULTA;
		
		/* modelo epidemiologico */
		
		GAMA		=	0.5*m_M*v[4]; //la mitad que madura son hembras
		
		b_theta_pV	=	bite_rate*theta_T(Tm)*MIObv;
		
		if ( Tdeath < NO_INFECCION ){ b_theta_pV = 0.;}
		
			/* parametros para la latencia del aedes*/
		
		mu_V		=	muerte_V(Tm)*MU_MOSQUITA_ADULTA;
		
		if( Tdeath < MATAR_VECTORES ){  mu_V = 2.*mu_V; }
			
		media_VE	=	0.1216*Tm*Tm - 8.66*Tm + 154.79;
				
		sigma_V		=	1./media_VE;
		
		
		EV 			=	sigma_V*v[7] - sigma_V*integral_1 + integral_2 - mu_V*integral_3 + integral_4;
		
		if ( EV < 0. ){ EV = 0.;}
		if(Tm<NO_LATENCIA){EV=0.;}
		
		b_theta_pH		=	bite_rate*theta_T(Tm)*MIObh;
		
		if ( Tdeath < NO_INFECCION ){ b_theta_pH = 0.;}
		
		sigma_H			=	1./Remove_expose;
		
		gama			=	1./Remove_infect;
		
		deltaI			=	importados[count];
		
		
		
		
		
		
		/* Resolucion de EDO */
		case_day 	=	 b_theta_pH*(v[9]/poblacion)*v[8];
		
		rk4(v,dv,m,t,dt,v,sistema);
		
		/* muerte a las fases acuaticas */
		if (v[0] < 0.){ v[0] = 0.;} //E_D
		
		if (v[1] < 0.){ v[1] = 0.;} //E_W
		if ( Tdeath < 10. ){ v[1] =  SUV*v[1];}
		
		if (v[2] < 0.){ v[2] = 0.;} //Larv
		if (v[2] > KL){ v[2] = KL;}
		if ( Tdeath < 10. ){ v[2] =  SUV*v[2];}
		
		if (v[3] < 0.){ v[3] = 0.;} //Pupa
		if ( Tdeath < 10. ){v[3] =  SUV*v[3] ;}
		
		if (v[4] < 0.){ v[4] = 0.;} //M
		
		if (v[5] < 0.){ v[5] = 0.;}
		
		if (v[6] < 0.){ v[6] = 0.;}
		
		if (v[7] < 0.){ v[7] = 0.;}
		
		if (v[8] < 0.){ v[8] = 0.;}
		
		if (v[9] < 0.){ v[9] = 0.;}
		
		if (v[10] < 0.){ v[10] = 0.;}
		
		if (v[11] < 0.){ v[11] = 0.;}
		
		if (v[12] < 0.){ v[12] = 0.;}
		//printf ("%d\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\n",count,KL, m_E_C_G, C_L, v[2]);
		
		//if (t > 40){break;}

		t+=dt;
		
		count++;
		
		day++;
		
		if ( day < 7){
			
			case_week 	= case_week + case_day;
			
			v_h_week   = v_h_week + b_theta_pH*b_theta_pV*( 1./(mu_V*gama ) )*pow( ( k_VE*sigma_V/(mu_V + k_VE*sigma_V) )  , k_VE)*( v[5]/poblacion );//*( v[9]/poblacion );      
			
			pop_Vect   = pop_Vect + v[5]/poblacion;
		}
		else {
			
			fprintf(casos,"%d\t%.2f\n",week,case_week);
			
			fprintf(vector_host,"%d\t%.2f\n",week, v_h_week/7.);
			
			fprintf(mosquitos,"%d\t%.2f\n",week, pop_Vect/7.);
			
			case_week	=	0.;
			
			v_h_week    =   0.;
			
			pop_Vect	= 	0.;
			
			day			=	0;
			week		=	week + 1;
		
			
			
		}
			
	/*	float nat;
		for(t=1.; t<=DAYS; t+=1)
		{
		nat= ((b_theta_pH * b_theta_pV)/(mu_V * gama))* (sigma_V/(mu_V+sigma_V)) * (v[5]/poblacion);
		fprintf(R_o,"%f\t%.2f\n",t,v[5]);
	}*/
	}
	
	
	
	
	return 0;
}

void sistema(double t, double v[], double dv[]){
	
	double E_D, E_W, L, P, M, V;
	
	double H, H_S, H_E, H_I, H_R, V_S, V_E, V_I;
		
	E_D 	=	v[0];
	E_W		=	v[1];
	L		=	v[2];
	P		=	v[3];
	M		=	v[4];
	V		=	v[5];
	V_S		=	v[6];
	V_E		=	v[7];
	V_I		=	v[8];
	H_S		=	v[9];
	H_E		=	v[10];
	H_I		=	v[11];
	H_R		=	v[12];
	H		=	H_S + H_E + H_I + H_R;
	
	
	dv[0]	=	beta_day_theta_0*V - fR*E_D - mu_Dry*E_D;
	
	dv[1]	=	fR*E_D - m_E_C_G*E_W - mu_Wet*E_W;
	
	dv[2]	=	m_E_C_G*E_W - m_L*L - ( mu_L + C_L )*L;
	
	dv[3]	=	m_L*L - m_P*P - mu_P*P;
	
	dv[4]	=	m_P*P - m_M*M - mu_M*M;
	
	dv[5]	=	GAMA - mu_V*V;
	
	dv[6]	=	GAMA - b_theta_pV*(H_I/H)*V_S - mu_V*V;
	
	dv[7]	=	b_theta_pV*(H_I/H)*V_S - EV - mu_V*V_E;
	
	dv[8]	=	EV - mu_V*V_I;
	
	dv[9]	=	- b_theta_pH*(H_S/H)*V_I - sigma_H*H_E;
	
	dv[10]	=	b_theta_pH*(H_S/H)*V_I - sigma_H*H_E;
	
	dv[11]	=	sigma_H*H_E - gama*H_I + deltaI;
	
	dv[12]	=	gama*H_I;
	
	return;
	
}

/*pdf gamma*/
double suv_exp(double tt,double s, double rate)
{
	double salida,
		tau = tt-s;
	
	if (tt<=0){ salida =0.;}
	
	else{
		salida =  exp(- rate * tau);
		}
	
	return salida;
}
/*Fbar - supervivencia*/
double Fbar(double tt,double s, double k, double theta)
{
	double salida,
		tau	=	tt-s;
	if (tt<0.)
	{
		salida=0.;
	}
	else
	{
		//salida = exp(-tau/(k*theta)); // distribucion exponencial
		
		salida = gammq(k,tau/theta); // distribucion gamma
		
		/*if (tau<k*theta){salida=1.;} // distribucion escalon
		else{salida=0.;}*/
	}
	return salida;
}
