#include "Parameter_files/INIT_PARAMS.H"
#include "Parameter_files/ANAL_PARAMS.H"


void HII_filter(fftwf_complex *box, int filter_type, float R){
  int n_x, n_z, n_y;
  float k_x, k_y, k_z, k_mag, k_mag_x, k_mag_y,kR;

#pragma omp parallel shared(box, filter_type, R) private(k_x, k_y, k_z, k_mag, k_mag_x, k_mag_y, kR, n_x, n_z, n_y)
{

  // loop through k-box
#pragma omp for
    for (n_x=0; n_x<HII_DIM; n_x++){
        if (n_x>HII_MIDDLE) {
//            k_x =(n_x-HII_DIM) * DELTA_K;
            k_mag_x = (n_x-HII_DIM)*(n_x-HII_DIM);
        }
        else {
//            k_x = n_x * DELTA_K;
            k_mag_x = n_x*n_x;
        }

//        k_mag_x = k_x*k_x;

        for (n_y=0; n_y<HII_DIM; n_y++){
            if (n_y>HII_MIDDLE) {
//                k_y =(n_y-HII_DIM) * DELTA_K;
                k_mag_y = (n_y-HII_DIM)*(n_y-HII_DIM);
            }
            else {
//                k_y = n_y * DELTA_K;
                k_mag_y = n_y*n_y;
            }

//            k_mag_y = k_y*k_y;

            for (n_z=0; n_z<=HII_MIDDLE; n_z++){
//                k_z = n_z * DELTA_K;

//                k_mag = sqrt(k_mag_x + k_mag_y + k_z*k_z);
                k_mag = sqrt(k_mag_x + k_mag_y + n_z*n_z)*DELTA_K;

                kR = k_mag*R; // real space top-hat
                if (filter_type == 0){ // real space top-hat
                    if (kR > 1e-4){
//                        box[HII_C_INDEX(n_x, n_y, n_z)] *= 3.0 * (sin(kR)/pow(kR, 3) - cos(kR)/pow(kR, 2));
                        box[HII_C_INDEX(n_x, n_y, n_z)] *= 3.0 * (sin(kR)/(kR*kR*kR) - cos(kR)/(kR*kR));
                    }
                }
                else if (filter_type == 1){ // k-space top hat
                    kR *= 0.413566994; // equates integrated volume to the real space top-hat (9pi/2)^(-1/3)
                    if (kR > 1){
                        box[HII_C_INDEX(n_x, n_y, n_z)] = 0;
                    }
                }
                else if (filter_type == 2){ // gaussian
                    kR *= 0.643; // equates integrated volume to the real space top-hat
                    box[HII_C_INDEX(n_x, n_y, n_z)] *= pow(E, -kR*kR/2.0);
                }
                else{
                    if ( (n_x==0) && (n_y==0) && (n_z==0) )
                        fprintf(stderr, "HII_filter.c: Warning, filter type %i is undefined\nBox is unfiltered\n", filter_type);
                }
            }
        }
    } // end looping through k box

}
 return;
}


/*
  all lengths are in units of the box size
  (x,y,z) is the closest reflection of (x2,y2,z2) to (x1, y1, z1)
*/
float distance_coord(float x1, float y1, float z1,
		     float x2, float y2, float z2,
		     float *x, float *y, float *z
		     ){
  float minimumsq, xsq, ysq, zsq, xplussq, yplussq, zplussq, xminsq, yminsq, zminsq;

    // remember to check all reflections
    xsq = pow(x1-x2, 2);
    ysq = pow(y1-y2, 2);
    zsq = pow(z1-z2, 2);
    xplussq = pow(x1-x2+1, 2);
    yplussq = pow(y1-y2+1, 2);
    zplussq = pow(z1-z2+1, 2);
    xminsq = pow(x1-x2-1, 2);
    yminsq = pow(y1-y2-1, 2);
    zminsq = pow(z1-z2-1, 2);

    // the true distance is the minimum of the permutations above
    minimumsq = 10;
    if ( minimumsq > (xsq + ysq + zsq) ){
      minimumsq = (xsq + ysq + zsq);
      *x=x2; *y=y2; *z=z2;
    }
    if (minimumsq > (xsq + ysq + zplussq)){
      minimumsq = (xsq + ysq + zplussq);
      *x=x2; *y=y2; *z=z2-1;
    }
    if (minimumsq > (xsq + ysq + zminsq)){
      minimumsq = (xsq + ysq + zminsq);
      *x=x2; *y=y2; *z=z2+1;
    }

    if (minimumsq > (xsq + yplussq + zsq)){
      minimumsq = (xsq + yplussq + zsq);
      *x=x2; *y=y2-1; *z=z2;
    }
    if (minimumsq > (xsq + yplussq + zplussq)){
      minimumsq = (xsq + yplussq + zplussq);
      *x=x2; *y=y2-1; *z=z2-1;
    }
    if (minimumsq > (xsq + yplussq + zminsq)){
      minimumsq = (xsq + yplussq + zminsq);
      *x=x2; *y=y2-1; *z=z2+1;
    }

    if (minimumsq > (xsq + yminsq + zsq)){
      minimumsq = (xsq + yminsq + zsq);
      *x=x2; *y=y2+1; *z=z2;
    }
    if (minimumsq > (xsq + yminsq + zplussq)){
      minimumsq = (xsq + yminsq + zplussq);
      *x=x2; *y=y2+1; *z=z2-1;
    }
    if (minimumsq > (xsq + yminsq + zminsq)){
      minimumsq = (xsq + yminsq + zminsq);
      *x=x2; *y=y2+1; *z=z2+1;
    }

    if (minimumsq > (xplussq + ysq + zsq)){
      minimumsq = (xplussq + ysq + zsq);
      *x=x2-1; *y=y2; *z=z2;
    }
    if (minimumsq > (xplussq + ysq + zplussq)){
      minimumsq = (xplussq + ysq + zplussq);
      *x=x2-1; *y=y2; *z=z2-1;
    }
    if (minimumsq > (xplussq + ysq + zminsq)){
      minimumsq = (xplussq + ysq + zminsq);
      *x=x2-1; *y=y2; *z=z2+1;
    }

    if (minimumsq > (xplussq + yplussq + zsq)){
      minimumsq = (xplussq + yplussq + zsq);
      *x=x2-1; *y=y2-1; *z=z2;
    }
    if (minimumsq > (xplussq + yplussq + zplussq)){
      minimumsq = (xplussq + yplussq + zplussq);
      *x=x2-1; *y=y2-1; *z=z2-1;
    }
    if (minimumsq > (xplussq + yplussq + zminsq)){
      minimumsq = (xplussq + yplussq + zminsq);
      *x=x2-1; *y=y2-1; *z=z2+1;
    }

    if (minimumsq > (xplussq + yminsq + zsq)){
      minimumsq = (xplussq + yminsq + zsq);
      *x=x2-1; *y=y2+1; *z=z2;
    }
    if (minimumsq > (xplussq + yminsq + zplussq)){
      minimumsq = (xplussq + yminsq + zplussq);
      *x=x2-1; *y=y2+1; *z=z2-1;
    }
    if (minimumsq > (xplussq + yminsq + zminsq)){
      minimumsq = (xplussq + yminsq + zminsq);
      *x=x2-1; *y=y2+1; *z=z2+1;
    }

    if (minimumsq > (xminsq + ysq + zsq)){
      minimumsq = (xminsq + ysq + zsq);
      *x=x2+1; *y=y2; *z=z2;
    }
    if (minimumsq > (xminsq + ysq + zplussq)){
      minimumsq = (xminsq + ysq + zplussq);
      *x=x2+1; *y=y2; *z=z2-1;
    }
    if (minimumsq > (xminsq + ysq + zminsq)){
      minimumsq = (xminsq + ysq + zminsq);
      *x=x2+1; *y=y2; *z=z2+1;
    }

    if (minimumsq > (xminsq + yplussq + zsq)){
      minimumsq = (xminsq + yplussq + zsq);
      *x=x2+1; *y=y2-1; *z=z2;
    }
    if (minimumsq > (xminsq + yplussq + zplussq)){
      minimumsq = (xminsq + yplussq + zplussq);
      *x=x2+1; *y=y2-1; *z=z2-1;
    }
    if (minimumsq > (xminsq + yplussq + zminsq)){
      minimumsq = (xminsq + yplussq + zminsq);
      *x=x2+1; *y=y2-1; *z=z2+1;
    }

    if (minimumsq > (xminsq + yminsq + zsq)){
      minimumsq = (xminsq + yminsq + zsq);
      *x=x2+1; *y=y2+1; *z=z2;
    }
    if (minimumsq > (xminsq + yminsq + zplussq)){
      minimumsq = (xminsq + yminsq + zplussq);
      *x=x2+1; *y=y2+1; *z=z2-1;
    }
    if (minimumsq > (xminsq + yminsq + zminsq)){
      minimumsq = (xminsq + yminsq + zminsq);
      *x=x2+1; *y=y2+1; *z=z2+1;
    }

    return sqrt(minimumsq);
}

/*
  all lengths are in units of the box size
*/
float distance(float x1, float y1, float z1, float x2, float y2, float z2){
  float minimumsq, xsq, ysq, zsq, xplussq, yplussq, zplussq, xminsq, yminsq, zminsq;

    // remember to check all reflections
    xsq = pow(x1-x2, 2);
    ysq = pow(y1-y2, 2);
    zsq = pow(z1-z2, 2);
    xplussq = pow(x1-x2+1, 2);
    yplussq = pow(y1-y2+1, 2);
    zplussq = pow(z1-z2+1, 2);
    xminsq = pow(x1-x2-1, 2);
    yminsq = pow(y1-y2-1, 2);
    zminsq = pow(z1-z2-1, 2);

    // the true distance is the minimum of the permutations above
    minimumsq = 10;
    if ( minimumsq > (xsq + ysq + zsq) )
      minimumsq = (xsq + ysq + zsq);
    if (minimumsq > (xsq + ysq + zplussq))
      minimumsq = (xsq + ysq + zplussq);
    if (minimumsq > (xsq + ysq + zminsq))
      minimumsq = (xsq + ysq + zminsq);

    if (minimumsq > (xsq + yplussq + zsq))
      minimumsq = (xsq + yplussq + zsq);
    if (minimumsq > (xsq + yplussq + zplussq))
      minimumsq = (xsq + yplussq + zplussq);
    if (minimumsq > (xsq + yplussq + zminsq))
      minimumsq = (xsq + yplussq + zminsq);

    if (minimumsq > (xsq + yminsq + zsq))
      minimumsq = (xsq + yminsq + zsq);
    if (minimumsq > (xsq + yminsq + zplussq))
      minimumsq = (xsq + yminsq + zplussq);
    if (minimumsq > (xsq + yminsq + zminsq))
      minimumsq = (xsq + yminsq + zminsq);


    if (minimumsq > (xplussq + ysq + zsq))
      minimumsq = (xplussq + ysq + zsq);
    if (minimumsq > (xplussq + ysq + zplussq))
      minimumsq = (xplussq + ysq + zplussq);
    if (minimumsq > (xplussq + ysq + zminsq))
      minimumsq = (xplussq + ysq + zminsq);

    if (minimumsq > (xplussq + yplussq + zsq))
      minimumsq = (xplussq + yplussq + zsq);
    if (minimumsq > (xplussq + yplussq + zplussq))
      minimumsq = (xplussq + yplussq + zplussq);
    if (minimumsq > (xplussq + yplussq + zminsq))
      minimumsq = (xplussq + yplussq + zminsq);

    if (minimumsq > (xplussq + yminsq + zsq))
      minimumsq = (xplussq + yminsq + zsq);
    if (minimumsq > (xplussq + yminsq + zplussq))
      minimumsq = (xplussq + yminsq + zplussq);
    if (minimumsq > (xplussq + yminsq + zminsq))
      minimumsq = (xplussq + yminsq + zminsq);


    if (minimumsq > (xminsq + ysq + zsq))
      minimumsq = (xminsq + ysq + zsq);
    if (minimumsq > (xminsq + ysq + zplussq))
      minimumsq = (xminsq + ysq + zplussq);
    if (minimumsq > (xminsq + ysq + zminsq))
      minimumsq = (xminsq + ysq + zminsq);

    if (minimumsq > (xminsq + yplussq + zsq))
      minimumsq = (xminsq + yplussq + zsq);
    if (minimumsq > (xminsq + yplussq + zplussq))
      minimumsq = (xminsq + yplussq + zplussq);
    if (minimumsq > (xminsq + yplussq + zminsq))
      minimumsq = (xminsq + yplussq + zminsq);

    if (minimumsq > (xminsq + yminsq + zsq))
      minimumsq = (xminsq + yminsq + zsq);
    if (minimumsq > (xminsq + yminsq + zplussq))
      minimumsq = (xminsq + yminsq + zplussq);
    if (minimumsq > (xminsq + yminsq + zminsq))
      minimumsq = (xminsq + yminsq + zminsq);

    return sqrt(minimumsq);
}


// helper function for update in sphere below
void check_region(float * box, int dimensions, float Rsq_curr_index, int x, int y, int z, int x_min, int x_max, int y_min, int y_max, int z_min, int z_max){
  int x_curr, y_curr, z_curr, x_index, y_index, z_index;
  float xsq, xplussq, xminsq, ysq, yplussq, yminsq, zsq, zplussq, zminsq;

  for (x_curr=x_min; x_curr<=x_max; x_curr++){
    for (y_curr=y_min; y_curr<=y_max; y_curr++){
      for (z_curr=z_min; z_curr<=z_max; z_curr++){
	x_index = x_curr;
	y_index = y_curr;
	z_index = z_curr;
	// adjust if we are outside of the box
	if (x_index<0) {x_index += dimensions;}
	else if (x_index>=dimensions) {x_index -= dimensions;}
	if (y_index<0) {y_index += dimensions;}
	else if (y_index>=dimensions) {y_index -= dimensions;}
	if (z_index<0) {z_index += dimensions;}
	else if (z_index>=dimensions) {z_index -= dimensions;}

	// now check
	//	printf("checking %i, %i, %i\n", x_index, y_index, z_index);
	//printf("this is index, %llu\n", HII_R_INDEX(x_index, y_index, z_index));
	//fflush(NULL);
	if (box[HII_R_INDEX(x_index, y_index, z_index)]){ // untaken pixel (not part of other halo)
	  //	  printf("in if\n");
	  // fflush(NULL);
	  // remember to check all reflections
	  xsq = pow(x-x_index, 2);
	  ysq = pow(y-y_index, 2);
	  zsq = pow(z-z_index, 2);
	  xplussq = pow(x-x_index+dimensions, 2);
	  yplussq = pow(y-y_index+dimensions, 2);
	  zplussq = pow(z-z_index+dimensions, 2);
	  xminsq = pow(x-x_index-dimensions, 2);
	  yminsq = pow(y-y_index-dimensions, 2);
	  zminsq = pow(z-z_index-dimensions, 2);
	  if ( (Rsq_curr_index > (xsq + ysq + zsq)) ||
	       (Rsq_curr_index > (xsq + ysq + zplussq)) ||
	       (Rsq_curr_index > (xsq + ysq + zminsq)) ||

	       (Rsq_curr_index > (xsq + yplussq + zsq)) ||
	       (Rsq_curr_index > (xsq + yplussq + zplussq)) ||
	       (Rsq_curr_index > (xsq + yplussq + zminsq)) ||

	       (Rsq_curr_index > (xsq + yminsq + zsq)) ||
	       (Rsq_curr_index > (xsq + yminsq + zplussq)) ||
	       (Rsq_curr_index > (xsq + yminsq + zminsq)) ||


	       (Rsq_curr_index > (xplussq + ysq + zsq)) ||
	       (Rsq_curr_index > (xplussq + ysq + zplussq)) ||
	       (Rsq_curr_index > (xplussq + ysq + zminsq)) ||

	       (Rsq_curr_index > (xplussq + yplussq + zsq)) ||
	       (Rsq_curr_index > (xplussq + yplussq + zplussq)) ||
	       (Rsq_curr_index > (xplussq + yplussq + zminsq)) ||

	       (Rsq_curr_index > (xplussq + yminsq + zsq)) ||
	       (Rsq_curr_index > (xplussq + yminsq + zplussq)) ||
	       (Rsq_curr_index > (xplussq + yminsq + zminsq)) ||


	       (Rsq_curr_index > (xminsq + ysq + zsq)) ||
	       (Rsq_curr_index > (xminsq + ysq + zplussq)) ||
	       (Rsq_curr_index > (xminsq + ysq + zminsq)) ||

	       (Rsq_curr_index > (xminsq + yplussq + zsq)) ||
	       (Rsq_curr_index > (xminsq + yplussq + zplussq)) ||
	       (Rsq_curr_index > (xminsq + yplussq + zminsq)) ||

	       (Rsq_curr_index > (xminsq + yminsq + zsq)) ||
	       (Rsq_curr_index > (xminsq + yminsq + zplussq)) ||
	       (Rsq_curr_index > (xminsq + yminsq + zminsq))
	     ){

	    // we are within the sphere defined by R, so change flag in box array
	    	    box[HII_R_INDEX(x_index, y_index, z_index)] = 0;
	    //	    box[HII_R_INDEX(x_index, y_index, z_index)] = 15;
	    //printf("%i, %i, %i\n", x_index, y_index, z_index);
		    //		    fflush(NULL);
	  }
	}
      }
    }
  }
}


/*
  Function UPDATE_IN_SPHERE takes in a <box> and flags all points
  which fall within radius R of (x,y,z).
  all lengths are in units of box size.
*/
void update_in_sphere(float * box, int dimensions, float R, float xf, float yf, float zf){
  int x_curr, y_curr, z_curr, xb_min, xb_max, yb_min, yb_max, zb_min, zb_max, R_index;
  int xl_min, xl_max, yl_min, yl_max, zl_min, zl_max;
  float Rsq_curr_index;
  int x_index, y_index, z_index, x, y, z;

  if (R<0) return;
  //  printf("in update, dim=%i, R=%f, x=%f, y=%f, z=%f\n", dimensions, R, xf, yf, zf);
  //fflush(NULL);
  // convert distances to index units
  x = (int) (xf * dimensions + 0.5); // +0.5 for rounding
  y = (int) (yf * dimensions + 0.5);
  z = (int) (zf * dimensions + 0.5);


  /*****  first, just automatically fill in the inner cube whose diagonal is R, side is 2R/sqrt(3) *****/
  R_index = ceil(R/sqrt(3.0)*dimensions)-1;
  // set parameter range
  xl_min = x-R_index;
  xl_max = x+R_index;
  yl_min = y-R_index;
  yl_max = y+R_index;
  zl_min = z-R_index;
  zl_max = z+R_index;

  for (x_curr=xl_min; x_curr<=xl_max; x_curr++){
    for (y_curr=yl_min; y_curr<=yl_max; y_curr++){
      for (z_curr=zl_min; z_curr<=zl_max; z_curr++){
	x_index = x_curr;
	y_index = y_curr;
	z_index = z_curr;
	// adjust if we are outside of the box
	if (x_index<0) {x_index += dimensions;}
	else if (x_index>=dimensions) {x_index -= dimensions;}
	if (y_index<0) {y_index += dimensions;}
	else if (y_index>=dimensions) {y_index -= dimensions;}
	if (z_index<0) {z_index += dimensions;}
	else if (z_index>=dimensions) {z_index -= dimensions;}

	// now just paint it
	//box[HII_R_INDEX(x_index, y_index, z_index)] = 15;
	       	box[HII_R_INDEX(x_index, y_index, z_index)] = 0;
      }
    }
  }


  /****** now check the pixels between the smaller and the larger cube which encloses the sphere ******/
  R_index = ceil(R*dimensions);
  Rsq_curr_index = pow(R*dimensions, 2); // convert to index
  // set parameter range
  xb_min = x-R_index;
  xb_max = x+R_index;
  yb_min = y-R_index;
  yb_max = y+R_index;
  zb_min = z-R_index;
  zb_max = z+R_index;

  //    check_region(box, dimensions, Rsq_curr_index, x,y,z, xb_min, xb_max, yb_min, yb_max, zb_min, zb_max);

  check_region(box, dimensions, Rsq_curr_index, x,y,z, xb_min, xl_min, yb_min, yb_max, zb_min, zb_max);
  check_region(box, dimensions, Rsq_curr_index, x,y,z, xl_max, xb_max, yb_min, yb_max, zb_min, zb_max);

  check_region(box, dimensions, Rsq_curr_index, x,y,z, xb_min, xb_max, yb_min, yl_min, zb_min, zb_max);
  check_region(box, dimensions, Rsq_curr_index, x,y,z, xb_min, xb_max, yl_max, yb_max, zb_min, zb_max);

  check_region(box, dimensions, Rsq_curr_index, x,y,z, xb_min, xb_max, yb_min, yb_max, zb_min, zl_min);
  check_region(box, dimensions, Rsq_curr_index, x,y,z, xb_min, xb_max, yb_min, yb_max, zl_max, zb_max);
}
