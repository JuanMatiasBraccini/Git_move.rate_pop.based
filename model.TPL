GLOBALS_SECTION
  #include <iostream>
  #include <fstream>
  #undef REPORT
  #define REPORT(object) report<<#object "\n" << object << endl;
  #include <time.h>
  time_t start,finish;
  long hour,minute,second;
  double elapsed_time;

DATA_SECTION
  init_int Tag_StartYear;  //first year of tagging
  init_int Tag_EndYear;    //last year of recaptures
  init_int Tag_nyrs;       //number of years
  init_int Tag_Areas;      //spatial strata. 1:West, 2:ZN1, 3:ZN2
  init_int Tag_numTagGp;   //number of tag groups
  init_3darray Tag_Releases_areas(1,Tag_numTagGp,1,1,1,3);  //numbers releases by tag group and region
  init_3darray Tag_Releases_yrs(1,Tag_numTagGp,1,1,1,1);  //release  years 
  init_3darray Tag_Recaptures(1,Tag_numTagGp,Tag_StartYear,Tag_EndYear,1,Tag_Areas);     //recaptures by tag group, (1st dim), years with recapture (2nd dim) and regions (3rd dim)
  init_3darray Tag_Recapture_index(1,Tag_numTagGp,Tag_StartYear,Tag_EndYear,1,Tag_Areas);     //index of recaptures 
  init_matrix Tag_Reporting(Tag_StartYear,Tag_EndYear,1,Tag_Areas);          // reporting rate by year and area
  init_number Tag_Shedding;                                      //tag shedding rate
  init_number Tag_M;                                             //natural mortality
  init_matrix Tag_Effort(Tag_StartYear,Tag_EndYear,1,Tag_Areas);              //Effort by year and area
  init_int Tag_nMovParams;                                        //number of movement parameters
  init_number liketype;                                       // select type of likelihood
  init_number likePoissonK;                                   // small constant
  //!!cout<<Tag_Recaptures<<endl;exit(1);
  
  int t; //declare indices
  int y;
  int a;
  int b;
  int Tag_yearTag;
  

    //!!cout<<"Tag_Releases_yrs\n"<<Tag_Releases_yrs(1)<<endl;exit(1);

PARAMETER_SECTION
  init_bounded_number movPar_1_1(0,1,1);          //movement parameters
  init_bounded_number movPar_1_2(0,1,2);
  init_bounded_number movPar_2_1(0,1,1);
  init_bounded_number movPar_2_2(0,1,1);
  init_bounded_number movPar_3_2(0,1,1);
  init_bounded_number movPar_3_3(0,1,2);
  
  init_number Tag_log_Q1(3);    //catchability by area
  init_number Tag_log_Q2(3);
  init_number Tag_log_Q3(3);

  init_bounded_number Tag_log_Tau(0.0001,10,1);    //overdispersion of negative binomial
  
  

  objective_function_value f

  //Declare objects used in Tag_Pop_dyn()
  matrix MovMat(1,Tag_Areas,1,Tag_Areas);		// Matrix of movement parameters 
  3darray Tag_N(1,Tag_numTagGp,Tag_StartYear,Tag_EndYear,1,Tag_Areas);  //Predicted numbers in populations for tag group, year and region
  3darray Tag_Pred_Rec(1,Tag_numTagGp,Tag_StartYear,Tag_EndYear,1,Tag_Areas);  //Predicted recaptures in populations for tag group, year and region
  matrix Tag_F(Tag_StartYear,Tag_EndYear,1,Tag_Areas);
  matrix Tag_Q(Tag_StartYear,Tag_EndYear,1,Tag_Areas);
  number Tag_Tau;
  3darray Tag_N_no_mov(1,Tag_numTagGp,Tag_StartYear,Tag_EndYear,1,Tag_Areas);  //dummy for storing N without movement



  //Declare objects used in Tag_observe()
  

  //Declare objects used in Tag_objective_fn()
  number Tag_NLL
  

PRELIMINARY_CALCS_SECTION

  

PROCEDURE_SECTION
  Tag_Pop_dyn();
  Tag_observe();
  Tag_objective_fn();
  

FUNCTION Tag_Pop_dyn
//note: compute survival and movement of each tag group

 // put catchability in normal space
 for(t=Tag_StartYear;t<=Tag_EndYear;t++)
 {
   Tag_Q(t,1)=mfexp(Tag_log_Q1);
   Tag_Q(t,2)=mfexp(Tag_log_Q2);
   Tag_Q(t,3)=mfexp(Tag_log_Q3);
 }


 // calculate fishing mortality at time and area
  for(t=Tag_StartYear;t<=Tag_EndYear;t++)
 {
   Tag_F(t)=elem_prod(Tag_Q(t),Tag_Effort(t));
 }

   //fill in Movement matrix
  MovMat(1,1)=movPar_1_1;
  MovMat(1,2)=movPar_1_2;
  MovMat(1,3)=1-(movPar_1_1+movPar_1_2);

  MovMat(2,1)=movPar_2_1;
  MovMat(2,2)=movPar_2_2;
  MovMat(2,3)=1-(movPar_2_1+movPar_2_2);

  MovMat(3,1)=1-(movPar_3_3+movPar_3_2);
  MovMat(3,2)=movPar_3_2;
  MovMat(3,3)=movPar_3_3;
  
 // fill in population numbers
 Tag_N.initialize();
 Tag_Pred_Rec.initialize();

 for(t=1; t<=Tag_numTagGp;t++)                                            //loop over tag groups
 {
    Tag_yearTag=Tag_Releases_yrs(t,1,1);
    for(y=Tag_yearTag; y<=Tag_EndYear;y++)                               //loop over years
    {         
       if(y==Tag_yearTag)
       {                  
         //Recruitment tags
         Tag_N(t,y)=Tag_Releases_areas(t,1);

         //Apply mortality
         for(a=1; a<=Tag_Areas;a++)                                    //loop over areas
          {
            Tag_N(t,y,a) =  Tag_N(t,y,a) * mfexp(-(Tag_F(y,a) + Tag_M + Tag_Shedding));
          }
          
        }else
        {                 
          //Apply mortality
          for(a=1; a<=Tag_Areas;a++)                                    //loop over areas
          {
           Tag_N(t,y,a) =  Tag_N(t,y-1,a) * mfexp(-(Tag_F(y,a) + Tag_M + Tag_Shedding));
          }   
        }

        //track population change when no movement
        Tag_N_no_mov(t,y)=Tag_N(t,y);
        
        //Apply movement
        Tag_N(t,y)=Tag_N(t,y)*MovMat;
    }
 }
             
             
 
FUNCTION Tag_observe
//note: 1. predict the expected number of tag recoveries using Baranov catch equation
//      2. compare to the actual number of tags that were captured and reported

 Tag_Pred_Rec.initialize();
 
 for(t=1; t<=Tag_numTagGp;t++)                               //loop over tag groups
 {
    Tag_yearTag=Tag_Releases_yrs(t,1,1);
    for(y=Tag_yearTag; y<=Tag_EndYear;y++)                   //loop over years
    {
      for(a=1; a<=Tag_Areas;a++)                            //loop over areas
      {
       Tag_Pred_Rec(t,y,a)=Tag_Reporting(y,a) * Tag_N(t,y,a) * ((Tag_F(y,a))/(Tag_F(y,a) + Tag_M + Tag_Shedding)) *(1-mfexp(-(Tag_F(y,a) + Tag_M + Tag_Shedding)));
      }
    }
 }

  // cout<<"Tag_N_no_mov\n"<<Tag_N_no_mov<<endl;
   //cout<<"Tag_N\n"<<Tag_N<<endl;
   //cout<<"Tag_Pred_Rec\n"<<Tag_Pred_Rec<<endl;exit(1);

FUNCTION Tag_objective_fn
//use a negative binomial error as it allows for substantial variability among observations
  Tag_Tau=mfexp(Tag_log_Tau);
      
  Tag_NLL=0;
  //cout<<"Tag_NLL "<<Tag_NLL<<endl;
  
 //Calculate observations' likelihood
 for(t=1; t<=Tag_numTagGp;t++)                               //loop over tag groups
 {
   for(y=Tag_StartYear; y<=Tag_EndYear;y++)                   //loop over years
   {
     for(a=1; a<=Tag_Areas;a++)                            //loop over areas
     {
       if(Tag_Recapture_index(t,y,a)>0)
       {
        // Least squares likelihood
 	if(liketype==1) Tag_NLL += square((Tag_Recaptures(t,y,a)+likePoissonK)-((Tag_Pred_Rec(t,y,a) + likePoissonK)));
	 
	// Poisson likelihood
	if(liketype==2) Tag_NLL += Tag_Pred_Rec(t,y,a) + likePoissonK - Tag_Recaptures(t,y,a)*log(Tag_Pred_Rec(t,y,a) + likePoissonK);
	 
        //Negative Binomial
        if(liketype==3) Tag_NLL +=(-log_negbinomial_density((Tag_Recaptures(t,y,a)+ likePoissonK),(Tag_Pred_Rec(t,y,a)+ likePoissonK),Tag_Tau));
       //cout<<"y "<<y<<" a "<<a <<"  Tag_NLL "<<Tag_NLL<<endl; 
       }
       
     }
   }
 }
 
 // Objective function	 
 f = Tag_NLL ;
 
 cout<<"f "<<f<<endl;





REPORT_SECTION
  REPORT(Tag_Pred_Rec);               //predicted recaptures


