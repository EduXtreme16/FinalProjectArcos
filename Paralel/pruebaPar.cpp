#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include <math.h>  
#include <vector>
#include <bits/stdc++.h> 
#include <chrono>
#include <omp.h>
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
	if (pendiente < -1.0 || pendiente > 1.0)
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

double Calculate_forcesx (double weighta, double weightb, double distancia){

	//cout << "Fuerza x con parametros: " << gravity <<" " << weighta << " " << weightb << " " << distancia << endl;

	double force = gravity*weighta*weightb/(distancia*distancia);
	if (force > 200)
		return 200;
	//cout << "Resultado: " << force << endl;
	return force;
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

double Calculate_forcesy (double weighta, double weightb, double distancia){
	
	double force = gravity*weighta*weightb/pow(distancia, 2) ;
	if (force > 200)
		return 200;
	return force;
}

double Calculate_acceleration ( vector <double> forces, double weight ){
	
	double acceleration = 0;

	for (unsigned int i = 0; i < forces.size(); i++){
		if (forces[i] != 0)
			acceleration += forces[i];
		//acceleration += forces.back();
		//forces.pop_back();
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
	
	using namespace std;
	using namespace std::chrono;
	using clk = chrono::high_resolution_clock;

	auto t1 = clk :: now();

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
		Asteroides[i].forcesx.resize(num_asteroides+num_planetas);
		Asteroides[i].forcesy.resize(num_asteroides+num_planetas);
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

	vector <int> CheckCrash (num_asteroides);


	double distancia, pendiente, angle, fuerza;
	//auto it = Asteroides[0].forcesx.end();

	omp_set_num_threads(8);



	/* Calculo de las fuerzas para cada par de elementos*/
	for (int times = 0; times < num_iteraciones ; times++){

		/* !!!! CALCULO DE FUERZAS DE ASTEROIDE i CON LOS DEMAS ASTEROIDES !!!! */

		for (unsigned int i = 0; i < Asteroides.size(); i++) {
			for (unsigned int j = 0; j < Asteroides[i].forcesx.size(); j++){
				Asteroides[i].forcesx[j] = 0;
			}
			for (unsigned int w = 0; w < Asteroides[i].forcesy.size(); w++){
				Asteroides[i].forcesy[w] = 0;			
			}
		}


		// inicializar los vectores de fuerzas a 0.

		fs3 << "--- asteroids vs asteroids ---" << endl;

		#pragma omp parallel for private(distancia, pendiente, angle, fuerza)
		for (unsigned int i = 0; i < Asteroides.size() ; i++){
			//vector<asteroid>::iterator it;
						
			for (unsigned int j = i+1; j < Asteroides.size(); j++){
				distancia = Calculate_dist(Asteroides[i].posx, Asteroides[i].posy, Asteroides[j].posx, Asteroides[j].posy);
				if (distancia > dmin){
					pendiente = Calculate_slope(Asteroides[i].posx, Asteroides[i].posy, Asteroides[j].posx, Asteroides[j].posy);
					angle = Calculate_angle(pendiente);


					//double fuerza = Calculate_forcesx(Asteroides[i].posx, Asteroides[i].posy, Asteroides[i].weight, Asteroides[j].posx, Asteroides[j].posy, Asteroides[j].weight );
					fuerza = Calculate_forcesx(Asteroides[i].weight, Asteroides[j].weight, distancia);
					
					if (i == 0 && j == 1)
						cout << "D:" << distancia << " " << "P:" << pendiente << " " << "A:" << angle << " " << "F:" << fuerza << endl;
					fs3 << i << " " << j << " " << fuerza << " " << angle << endl;
					/*omp_set_lock(&l);
					it = Asteroides[i].forcesx.end();
					Asteroides[i].forcesx.insert(it, fuerza*cos(angle));
					it = Asteroides[j].forcesx.end();
					Asteroides[j].forcesx.insert(it, -fuerza*cos(angle));
					*/

					Asteroides[i].forcesx[j] = fuerza*cos(angle);
					Asteroides[j].forcesx[i] = -fuerza*cos(angle);
					//fuerza = Calculate_forcesy(Asteroides[i].posx, Asteroides[i].posy, Asteroides[i].weight, Asteroides[j].posx, Asteroides[j].posy, Asteroides[j].weight );

					//fuerza = Calculate_forcesy(Asteroides[i].weight, Asteroides[j].weight, distancia);
					/*it = Asteroides[i].forcesy.end();
					Asteroides[i].forcesy.insert(it, fuerza*sin(angle));
					it = Asteroides[j].forcesy.end();
					Asteroides[j].forcesy.insert(it, -fuerza*sin(angle));
					omp_unset_lock(&l);*/
					Asteroides[i].forcesy[j] = fuerza*sin(angle);
					Asteroides[j].forcesy[i] = -fuerza*sin(angle);
				}
			}
		}

		fs3 << "--- asteroids vs planets ---" << endl;
		
		#pragma omp parallel for private(distancia, pendiente, angle, fuerza)
		/* !!!! CALCULO DE FUERZAS DE ASTEROIDE i CON LOS PLANETAS !!!! */
		for (unsigned int i = 0; i < Asteroides.size() ; i++){
			for (unsigned int j = 0; j < Planetas.size(); j++){
				distancia = Calculate_dist(Asteroides[i].posx, Asteroides[i].posy, Planetas[j].posx, Planetas[j].posy);
				//if (distancia > 2){

					pendiente = Calculate_slope(Asteroides[i].posx, Asteroides[i].posy, Planetas[j].posx, Planetas[j].posy);
					angle = Calculate_angle(pendiente);

					fuerza = Calculate_forcesx(Asteroides[i].weight, Planetas[j].weight, distancia);
					//double fuerza = Calculate_forcesx(Asteroides[i].posx, Asteroides[i].posy, Asteroides[i].weight, Planetas[j].posx, Planetas[j].posy, Planetas[j].weight );

					fs3 << i << " " << j << " " << fuerza << " " << angle  << endl;
					/*it = Asteroides[i].forcesx.end();
					Asteroides[i].forcesx.insert(it, fuerza*cos(angle));*/

					Asteroides[i].forcesx[num_asteroides+j] = fuerza*cos(angle);

					//fuerza = Calculate_forcesy(Asteroides[i].weight, Planetas[j].weight, distancia);
					//fuerza = Calculate_forcesy(Asteroides[i].posx, Asteroides[i].posy, Asteroides[i].weight, Planetas[j].posx, Planetas[j].posy, Planetas[j].weight );

					/*it = Asteroides[i].forcesy.end();
					Asteroides[i].forcesy.insert(it, fuerza*sin(angle));*/

					Asteroides[i].forcesy[num_asteroides+j] = fuerza*sin(angle);

				//}
			}
		}
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

		/* !!!! CALCULO DEL MOVIMIENTO !!!! */

		for (unsigned int i = 0; i < Asteroides.size() ; i++){
			/* Cálculo de los distintos parámetros que necesitamos (aceleración, velocidad y la posición final) */
			
			Asteroides[i].accelerationx = Calculate_acceleration(Asteroides[i].forcesx, Asteroides[i].weight);
			Asteroides[i].speedx = Calculate_speed(Asteroides[i].speedx, Asteroides[i].accelerationx);
			Asteroides[i].posx = Calculate_position(Asteroides[i].posx, Asteroides[i].speedx);


			Asteroides[i].accelerationy = Calculate_acceleration(Asteroides[i].forcesy, Asteroides[i].weight);
			Asteroides[i].speedy = Calculate_speed(Asteroides[i].speedy, Asteroides[i].accelerationy);
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

		for (unsigned int i = 0; i < CheckCrash.size(); i++){
			CheckCrash[i] = 0;
		}

		/* !!!! REBOTE CON ASTEROIDES !!!! */
		for (unsigned int i = 0; i < Asteroides.size(); i++){
			/*if (CheckCrash[i] == 0){
				CheckCrash[i] = 1;*/
				unsigned int aux = i;
				for (unsigned int j = i+1; j < Asteroides.size(); j++){
					distancia = Calculate_dist(Asteroides[i].posx, Asteroides[i].posy, Asteroides[j].posx, Asteroides[j].posy);
					if (distancia <= dmin && CheckCrash[j] == 0){
						//cout << "Chocan " << i <<" " << j << endl;
						double aux_speed = Asteroides[aux].speedx;
						Asteroides[aux].speedx = Asteroides[j].speedx;
						Asteroides[j].speedx = aux_speed;

						aux_speed = Asteroides[aux].speedy;
						Asteroides[aux].speedy = Asteroides[j].speedy;
						Asteroides[j].speedy = aux_speed;
						
						//CheckCrash[j] = 1;
						//aux = j;

					}
				}

			//}
		}

		fs3 << "\n ******************** ITERATION *******************" << endl;
	}

	fs3.close();

	ofstream fs2("out.txt"); 

	for (unsigned int i = 0; i < Asteroides.size() ; i++){
		fs2 << fixed << setprecision(3) << Asteroides[i].posx << " " << Asteroides[i].posy << " " << Asteroides[i].speedx << " " << Asteroides[i].speedy << " " << Asteroides[i].weight << endl;
	}

	fs2.close();
	auto t2 = clk::now();
	auto diff = duration_cast<microseconds>(t2-t1);
	cout << "HA TARDADO : " << diff.count() << "MICROSEGUNDOS" << endl;
	return 0;
}
