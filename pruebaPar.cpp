#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include <math.h>  
#include <vector>
#include <bits/stdc++.h>
#include <omp.h> 
#include <chrono>

using namespace std;

double gravity = 6.674e-5;
double timelapse = 0.1;
double dmin = 2.0;
double width = 200;
double height = 200;
double mass = 1000;
double sdm = 50;

struct asteroid{
	double posx;
	double posy;
	double weight;
	vector <double> forcesx;
	vector <double> forcesy;
	double speedx = 0;
	double speedy = 0;
	double accelerationx = 0;
	double accelerationy = 0;
} ;

struct planet{
	double posx;
	double posy;
	double weight;
} ;

double Calculate_dist (double xa, double ya, double xb, double yb){

	return sqrt(pow((xa-xb),2)+pow((ya-yb),2));
}

double Calculate_slope (double xa, double ya, double xb, double yb){

	double pendiente = (ya-yb)/(xa-xb);
	if (pendiente < -1 || pendiente > 1)
		pendiente = pendiente - (int(pendiente)/1);
	return pendiente;
}

double Calculate_angle(double xa, double ya, double xb, double yb ){

	return atan(Calculate_slope(xa, ya, xb, yb));
}

double Calculate_angle(double slope){
	return atan(slope);
}

double Calculate_forcesx (double xa, double ya, double weighta, double xb, double yb, double weightb){

	double distancia = Calculate_dist(xa, ya, xb, yb);

	double pendiente = Calculate_slope(xa, ya, xb, yb);
	double angle = Calculate_angle(pendiente);
	

	double force = gravity*weighta*weightb/pow(distancia, 2) * cos(angle);
	if (force > 200)
		force = 200;
	return force;
}

double Calculate_forcesx (double weighta, double weightb, double distancia, double angulo){

	cout << "Fuerza x con parametros: " << gravity <<" " << weighta << " " << weightb << " " << distancia << " " << angulo << endl;

	double force = gravity*weighta*weightb/(distancia*distancia);
	if (force > 200)
		return 200;
	cout << "Resultado: " << force << endl;
	return force * cos(angulo);
}

double Calculate_forcesy (double xa, double ya, double weighta, double xb, double yb, double weightb){

	//cout << "############### CÁLCULO DE LAS FUERZAS ############" << endl;
	//cout << "Parámetros: Gravity = " << gravity << "; weighta= " << weighta << "; weightb= " << weightb << "; Distancias al cuadrado= " << pow(Calculate_dist(xa, ya, xb, yb), 2) << "; Seno del angulo= " << sin(Calculate_angle(xa, ya, xb, yb)) << endl;
	double distancia = Calculate_dist(xa, ya, xb, yb);

	double pendiente = Calculate_slope(xa, ya, xb, yb);
	double angle = Calculate_angle(pendiente);

	double force = ((gravity*weighta*weightb)/(pow(distancia, 2))) * sin(angle);
	if (force > 200)
		force = 200;
	return force;
}

double Calculate_forcesy (double weighta, double weightb, double distancia, double angulo){
	
	double force = gravity*weighta*weightb/pow(distancia, 2) ;
	if (force > 200)
		return 200;
	return force * sin(angulo);
}

double Calculate_acceleration ( vector <double> forces, double weight ){
	
	double acceleration = 0;

	for (unsigned int i = 0; i < forces.size(); i++){
		acceleration += forces[i];
	}
	return acceleration / weight;
}

double Calculate_speed ( double speed, double acceleration){
	return speed + acceleration * timelapse;
}

double Calculate_position ( double position, double speed ){
	return position + speed * timelapse;
}

void printError(){
	cout << "nasteroids-seq: Wrong arguments.\nCorrect use:\nnasteroids-seq num_asteroides num_iteraciones num_planetas semilla" << endl;
}

int main(int argc, char** argv) {
	int pendiente = 0;
	int angle = 0;
	int fuerza = 0;
	omp_lock_t m;
	omp_init_lock(&m);
	// Cogemos el tiempo al inicio
	using namespace std::chrono;
	using clk = chrono::high_resolution_clock;
	auto t1 = clk::now();

	if (argc != 5){
		printError();
		return -1;
	}

	//cout << int(argv[1][0]) << endl;
	//cout << strlen(argv[1]) << endl;
	//cout << typeid(strlen(argv[1])).name() << endl;


	/* !!!! COMPROBACION DE PARAMETROS DE ENTRADA !!!! */
	for( unsigned int i = 1; i < 5; i++){
		if (atoi(argv[i]) > 2147483647 || atoi(argv[i]) < 0){
			printError();
			return -1;
		}
		for (unsigned int j = 0; j < strlen(argv[i]); j++){
			if (int(argv[i][j]) < 48 || int(argv[i][j]) > 57){
				printError();
				return -1;
			}
		}
	}

	int num_asteroides = atoi(argv[1]);
	int num_iteraciones = atoi(argv[2]);
	int num_planetas = atoi(argv[3]);
	int semilla = atoi(argv[4]);


	/* !!!! INICIALIZACION DEL ENGINE !!!! */
	default_random_engine re{semilla};
	uniform_real_distribution<double> xdist{0.0, std::nextafter(width,std :: numeric_limits<double>::max())};
	uniform_real_distribution<double> ydist{0.0, std::nextafter(height,std :: numeric_limits<double>::max())};
	normal_distribution<double> mdist{mass, sdm};

	//cout << num_asteroides << " " << num_iteraciones << " " << num_planetas << " " << semilla << endl; 

	/*ofstream fs("init_conf.txt"); 
	fs << num_asteroides << " " << num_iteraciones << " " << num_planetas << " " << semilla << endl;
	fs.close();
	*/

	/* !!!! CREACION DE LOS VECTORES DE ASTEROIDES Y PLANETAS !!!! */
	//asteroid Asteroides [num_asteroides];
	vector <asteroid> Asteroides (num_asteroides);
	

	//planet Planetas [num_planetas];
	vector <planet> Planetas (num_planetas);

	ofstream fs("init_conf.txt"); 
	fs << num_asteroides << " " << num_iteraciones << " " << num_planetas << " " << semilla << endl;


	//cout << "#### ATENSION KHE BIENEN LOS ASTIROIDES XDXD ####" << endl;

	/* !!!! COLOCACION DE ASTEROIDES, CON SU POSICION Y MASA !!!! */
	for (int i = 0; i < num_asteroides; i++){
		Asteroides[i].posx = xdist(re);
		Asteroides[i].posy = ydist(re);
		Asteroides[i].weight = mdist(re);
		//cout << Asteroides[i].posx << " " << Asteroides[i].posy << " " << Asteroides[i].weight << endl;
		fs << fixed << setprecision(3) << Asteroides[i].posx << " " << Asteroides[i].posy << " " << Asteroides[i].weight << endl;
	}

	//cout << "#### ATENSION KHE BIENEN LOS PLANETAS XDXD ####" << endl;

	/* !!!! CCOLOCACION DE LOS PLANETAS, CON SU POSICION Y MASA !!!! */
	for (int i = 0; i < num_planetas; i++){
		switch(i%4){
			case 0:
				Planetas[i].posx = 0;
				Planetas[i].posy = ydist(re);
				Planetas[i].weight = mdist(re)*10;
				break;
			case 1:
				Planetas[i].posx = xdist(re);
				Planetas[i].posy = 0;
				Planetas[i].weight = mdist(re)*10;
				break;

			case 2:
				Planetas[i].posx = width;
				Planetas[i].posy = ydist(re);
				Planetas[i].weight = mdist(re)*10;
				break;
			case 3:
				Planetas[i].posx = xdist(re);
				Planetas[i].posy = height;
				Planetas[i].weight = mdist(re)*10;
				break;
		}
		
		fs << Planetas[i].posx << " " << Planetas[i].posy << " " << Planetas[i].weight << endl;
		//cout << Planetas[i].posx << " " << Planetas[i].posy << " " << Planetas[i].weight << endl;
	}

	fs.close();

	ofstream fs3("step_by_step.txt"); 

	/* Calculo de las fuerzas para cada par de elementos*/
	for (int times = 0; times < num_iteraciones ; times++){
		fs3 << "******************** ITERATION *******************" << endl;

		/* !!!! CALCULO DE FUERZAS DE ASTEROIDE i CON LOS DEMAS ASTEROIDES !!!! */




		// inicializar los vectores de fuerzas a 0.

		fs3 << "--- asteroids vs asteroids ---" << endl;
		#pragma omp parallel num_threads(16)
		{
		int n = omp_get_num_threads();
		cout << "Número de Hilos : " <<  n << endl; 
		for (unsigned int i = 0; i < Asteroides.size() ; i++){
			//vector<asteroid>::iterator it;
			// Vamos a paralelizar el siguiente bucle con OpenMP, ya que es un bucle con demasiada cara
			#pragma omp for private(pendiente, angle, fuerza)
			for (unsigned int j = i+1; j < Asteroides.size(); j++){
				double distancia = Calculate_dist(Asteroides[i].posx, Asteroides[i].posy, Asteroides[j].posx, Asteroides[j].posy);
				if (distancia > 2){

					omp_set_lock(&m);
					pendiente = Calculate_slope(Asteroides[i].posx, Asteroides[i].posy, Asteroides[j].posx, Asteroides[j].posy);
					angle = Calculate_angle(pendiente);

					//omp_set_lock(&m);
					//double fuerza = Calculate_forcesx(Asteroides[i].posx, Asteroides[i].posy, Asteroides[i].weight, Asteroides[j].posx, Asteroides[j].posy, Asteroides[j].weight );
					fuerza = Calculate_forcesx(Asteroides[i].weight, Asteroides[j].weight, distancia, angle);

					fs3 << i << " " << j << " " << fuerza << " " << angle << endl;
					omp_unset_lock(&m);

					omp_set_lock(&m);
					auto it = Asteroides[i].forcesx.end();
					omp_unset_lock(&m);
	
					omp_set_lock(&m);
					Asteroides[i].forcesx.insert(it, fuerza);
					it = Asteroides[j].forcesx.end();
					Asteroides[j].forcesx.insert(it, -fuerza);
					//omp_unset_lock(&m);
					//fuerza = Calculate_forcesy(Asteroides[i].posx, Asteroides[i].posy, Asteroides[i].weight, Asteroides[j].posx, Asteroides[j].posy, Asteroides[j].weight );
					//omp_set_lock(&m);
					fuerza = Calculate_forcesy(Asteroides[i].weight, Asteroides[j].weight, distancia, angle);
					it = Asteroides[i].forcesy.end();
					Asteroides[i].forcesy.insert(it, fuerza);
					it = Asteroides[j].forcesy.end();
					Asteroides[j].forcesy.insert(it, -fuerza);
					omp_unset_lock(&m);
				}
			}
		}

		fs3 << "--- asteroids vs planets ---" << endl;

		/* !!!! CALCULO DE FUERZAS DE ASTEROIDE i CON LOS PLANETAS !!!! */
		#pragma omp for private(pendiente, angle, fuerza)

		for (unsigned int i = 0; i < Asteroides.size() ; i++){
//			#pragma omp for -> Da error , tengo que especificar el hilo maestro 
			for (unsigned int j = 0; j < Planetas.size(); j++){
				double distancia = Calculate_dist(Asteroides[i].posx, Asteroides[i].posy, Planetas[j].posx, Planetas[j].posy);
				//if (distancia > 2){
					
					
					omp_set_lock(&m);
					double pendiente = Calculate_slope(Asteroides[i].posx, Asteroides[i].posy, Planetas[j].posx, Planetas[j].posy);
					double angle = Calculate_angle(pendiente);

					double fuerza = Calculate_forcesx(Asteroides[i].weight, Planetas[j].weight, distancia, angle);
					//double fuerza = Calculate_forcesx(Asteroides[i].posx, Asteroides[i].posy, Asteroides[i].weight, Planetas[j].posx, Planetas[j].posy, Planetas[j].weight );

					fs3 << i << " " << j << " " << fuerza << " " << angle  << endl;
					auto it = Asteroides[i].forcesx.end();
					Asteroides[i].forcesx.insert(it, fuerza);

					fuerza = Calculate_forcesy(Asteroides[i].weight, Planetas[j].weight, distancia, angle);
					//fuerza = Calculate_forcesy(Asteroides[i].posx, Asteroides[i].posy, Asteroides[i].weight, Planetas[j].posx, Planetas[j].posy, Planetas[j].weight );

					it = Asteroides[i].forcesy.end();
					Asteroides[i].forcesy.insert(it, fuerza);
					omp_unset_lock(&m);

			}
		}
		omp_destroy_lock(&m);
			//cout << "#### Fuerzas de los asteroides en eje x ####" << endl;
			//for (auto x: Asteroides[i].forcesx) {
	        	//cout << x << "\t";
	    	//}
	    	//cout << endl;

	    	//cout << "#### Fuerzas de los asteroides en eje y ####" << endl;
			//for (auto x: Asteroides[i].forcesy) {
	        	//cout << x << "\t";
	    	//}
			//cout << endl;

		} // Hasta aquí vamos a paralelizar fin del #pragma omp parallel

		/* !!!! CALCULO DEL MOVIMIENTO !!!! */

		for (unsigned int i = 0; i < Asteroides.size() ; i++){
			/* Cálculo de los distintos parámetros que necesitamos (aceleración, velocidad y la posición final) */
			
			Asteroides[i].accelerationx = Calculate_acceleration(Asteroides[i].forcesx, Asteroides[i].weight);
			Asteroides[i].accelerationy = Calculate_acceleration(Asteroides[i].forcesy, Asteroides[i].weight);

			Asteroides[i].speedx = Calculate_speed(Asteroides[i].speedx, Asteroides[i].accelerationx);
			Asteroides[i].speedy = Calculate_speed(Asteroides[i].speedy, Asteroides[i].accelerationy);
			
			Asteroides[i].posx = Calculate_position(Asteroides[i].posx, Asteroides[i].speedx);
			Asteroides[i].posy = Calculate_position(Asteroides[i].posy, Asteroides[i].speedy);

			//cout << "Posicion final de Asteroide numero " << i << ": X:" << Asteroides[i].posx << " Y:" << Asteroides[i].posy << endl;
		}

		/* !!!! REBOTE CON BORDES !!!! */
		for (unsigned int i = 0; i < Asteroides.size() ; i++){
			if ((Asteroides[i].posx) >= height){
				Asteroides[i].posx = height -2;
				Asteroides[i].speedx = -Asteroides[i].speedx;
			}
			if ((Asteroides[i].posx) <= 0){
				Asteroides[i].posx = 2;
				Asteroides[i].speedx = -Asteroides[i].speedx;
			}
			if ((Asteroides[i].posy) >= width){
				Asteroides[i].posy = width-2;
				Asteroides[i].speedy = -Asteroides[i].speedy;
			}
			if ((Asteroides[i].posy) <= 0){
				Asteroides[i].posy = 2;
				Asteroides[i].speedy = -Asteroides[i].speedy;
			}
		}

		/* !!!! REBOTE CON ASTEROIDES !!!! */
		for (unsigned int i = 0; i < Asteroides.size() ; i++){
			unsigned int aux = i;
			for (unsigned int j = i+1; j < Asteroides.size(); j++){
				if (Calculate_dist(Asteroides[i].posx, Asteroides[i].posy, Asteroides[j].posx, Asteroides[j].posy) <= 2){
					cout << "Chocan " << i <<" " << j << endl;
					double aux_speed = Asteroides[aux].speedx;
					Asteroides[aux].speedx = Asteroides[j].speedx;
					Asteroides[j].speedx = aux_speed;

					aux_speed = Asteroides[aux].speedy;
					Asteroides[aux].speedy = Asteroides[j].speedy;
					Asteroides[j].speedy = aux_speed;

					aux = j;
				}
			}
		}
	}

	fs3.close();

	ofstream fs2("outP.txt"); 

	for (unsigned int i = 0; i < Asteroides.size() ; i++){
		fs2 << fixed << setprecision(3) << Asteroides[i].posx << " " << Asteroides[i].posy << " " << Asteroides[i].speedx << " " << Asteroides[i].speedy << " " << Asteroides[i].weight << endl;
	}

	fs2.close();

	// Una vez el programa haya acabado, tomamos este momento como referencia para
	// calcular lo que ha tardado.

	auto t2 = clk::now();
	auto diff = duration_cast<microseconds>(t2-t1);
	cout << "HA TARDADO : " << diff.count() << "MICROSEGUNDOS" << endl;

	return 0;
}
