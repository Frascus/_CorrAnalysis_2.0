#include "../include/virtual_FF_t_interval_list.h"


using namespace std;



void Get_virtual_ff_fit_interval(string W, double off, double tht, int &Tmin, int &Tmax) {


  if( W == "H1") {
    if( (tht > 0.2297 -eps(3)) && (tht < 0.2297 + eps(3)) && (off > 0.2003 - eps(3)) && (off < 0.2003 + eps(3))) { Tmin= 11; Tmax=17;     }
    else if( (tht > 0.2756 -eps(3)) && (tht < 0.2756 + eps(3)) && (off > 0.1687 - eps(3)) && (off < 0.1687 + eps(3))) {  Tmin= 13; Tmax=17;     }
    else if( (tht > 0.3603 -eps(3)) && (tht < 0.3603 + eps(3)) && (off > 0.1687 - eps(3)) && (off < 0.1687 + eps(3))) {  Tmin= 13; Tmax=17;     }
    else if( (tht > 0.2975 -eps(3)) && (tht < 0.2975 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 13; Tmax=18;     }
     else if( (tht > 0.4064 -eps(3)) && (tht < 0.4064 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 13; Tmax=18;     }
    else if( (tht > 0.4620 -eps(3)) && (tht < 0.4620 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 13; Tmax=18;     }
    else if( (tht > 0.3005 -eps(3)) && (tht < 0.3005 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 13; Tmax=20;     } //if mixed [8:20] //nj14 13-14  21-20
    else if( (tht > 0.4234 -eps(3)) && (tht < 0.4234 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 14; Tmax=20;     } //nj14 13-20
    else if( (tht > 0.4986 -eps(3)) && (tht < 0.4986 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 13; Tmax=20;     }
    else if( (tht > 0.5409 -eps(3)) && (tht < 0.5409 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 13; Tmax=20;     } //nj14 12-13 20-19
    else if( (tht > 0.2852 -eps(3)) && (tht < 0.2852 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 13; Tmax=20;     }
    else if( (tht > 0.4150 -eps(3)) && (tht < 0.4150 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 13; Tmax=20;     }
    else if( (tht > 0.5037 -eps(3)) && (tht < 0.5037 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 13; Tmax=20;     }
    else if( (tht > 0.5636 -eps(3)) && (tht < 0.5636 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 13; Tmax=20;     }
    else if( (tht > 0.5986 -eps(3)) && (tht < 0.5986 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 13; Tmax=19;     }
    else crash("In Get_virtual_ff_fit_interval I cannot find the kinematical point: off: "+to_string(off)+" , tht: "+to_string(tht));
  }

  else if (W == "H2") {
    if( (tht > 0.2297 -eps(3)) && (tht < 0.2297 + eps(3)) && (off > 0.2003 - eps(3)) && (off < 0.2003 + eps(3))) { Tmin= 8; Tmax=20;     }
    else if( (tht > 0.2756 -eps(3)) && (tht < 0.2756 + eps(3)) && (off > 0.1687 - eps(3)) && (off < 0.1687 + eps(3))) {  Tmin= 14; Tmax=20;     }
    else if( (tht > 0.3603 -eps(3)) && (tht < 0.3603 + eps(3)) && (off > 0.1687 - eps(3)) && (off < 0.1687 + eps(3))) {  Tmin= 14; Tmax=20;     }
    else if( (tht > 0.2975 -eps(3)) && (tht < 0.2975 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 15; Tmax=21;     }
    else if( (tht > 0.4064 -eps(3)) && (tht < 0.4064 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 15; Tmax=21;     } //nj24 15-20
    else if( (tht > 0.4620 -eps(3)) && (tht < 0.4620 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 15; Tmax=20;     }
    else if( (tht > 0.3005 -eps(3)) && (tht < 0.3005 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 18; Tmax=25;     } //nj24 17-25
    else if( (tht > 0.4234 -eps(3)) && (tht < 0.4234 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 15; Tmax=24;     }
    else if( (tht > 0.4986 -eps(3)) && (tht < 0.4986 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 14; Tmax=24;     }
    else if( (tht > 0.5409 -eps(3)) && (tht < 0.5409 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 14; Tmax=23;     }  //nj24 14-22
    else if( (tht > 0.2852 -eps(3)) && (tht < 0.2852 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 19; Tmax=25;     }
    else if( (tht > 0.4150 -eps(3)) && (tht < 0.4150 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 17; Tmax=24;     }
    else if( (tht > 0.5037 -eps(3)) && (tht < 0.5037 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 15; Tmax=24;     } //nj24 15-23
    else if( (tht > 0.5636 -eps(3)) && (tht < 0.5636 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 15; Tmax=23;     } //nj24 15-22
    else if( (tht > 0.5986 -eps(3)) && (tht < 0.5986 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 14; Tmax=22;     }
    else crash("In Get_virtual_ff_fit_interval I cannot find the kinematical point: off: "+to_string(off)+" , tht: "+to_string(tht));


  }
  
  else if (W == "FA") {
    if( (tht > 0.2297 -eps(3)) && (tht < 0.2297 + eps(3)) && (off > 0.2003 - eps(3)) && (off < 0.2003 + eps(3))) { Tmin= 13; Tmax=19;     }
    else if( (tht > 0.2756 -eps(3)) && (tht < 0.2756 + eps(3)) && (off > 0.1687 - eps(3)) && (off < 0.1687 + eps(3))) {  Tmin= 14; Tmax=19;     }
    else if( (tht > 0.3603 -eps(3)) && (tht < 0.3603 + eps(3)) && (off > 0.1687 - eps(3)) && (off < 0.1687 + eps(3))) {  Tmin= 14; Tmax=17;     }
    else if( (tht > 0.2975 -eps(3)) && (tht < 0.2975 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 14; Tmax=19;     }
    else if( (tht > 0.4064 -eps(3)) && (tht < 0.4064 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 14; Tmax=18;     } //doubt!
    else if( (tht > 0.4620 -eps(3)) && (tht < 0.4620 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 14; Tmax=18;     } //doubt!
    else if( (tht > 0.3005 -eps(3)) && (tht < 0.3005 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 6; Tmax=14;     }
    else if( (tht > 0.4234 -eps(3)) && (tht < 0.4234 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 8; Tmax=13;     }
    else if( (tht > 0.4986 -eps(3)) && (tht < 0.4986 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 9; Tmax=13;     }
    else if( (tht > 0.5409 -eps(3)) && (tht < 0.5409 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 9; Tmax=14;     }
    else if( (tht > 0.2852 -eps(3)) && (tht < 0.2852 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) { Tmin=10;Tmax=15;}  //Tmin= 10; Tmax=15;     }
    else if( (tht > 0.4150 -eps(3)) && (tht < 0.4150 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 7; Tmax=15;     }
    else if( (tht > 0.5037 -eps(3)) && (tht < 0.5037 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 8; Tmax=15;     }
    else if( (tht > 0.5636 -eps(3)) && (tht < 0.5636 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 9; Tmax=15;     }
    else if( (tht > 0.5986 -eps(3)) && (tht < 0.5986 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 9; Tmax=14;     }
    else crash("In Get_virtual_ff_fit_interval I cannot find the kinematical point: off: "+to_string(off)+" , tht: "+to_string(tht));



   }


  else if (W =="FV") {
    if( (tht > 0.2297 -eps(3)) && (tht < 0.2297 + eps(3)) && (off > 0.2003 - eps(3)) && (off < 0.2003 + eps(3))) { Tmin= 11; Tmax=19;     }
    else if( (tht > 0.2756 -eps(3)) && (tht < 0.2756 + eps(3)) && (off > 0.1687 - eps(3)) && (off < 0.1687 + eps(3))) {  Tmin= 11; Tmax=19;     }
    else if( (tht > 0.3603 -eps(3)) && (tht < 0.3603 + eps(3)) && (off > 0.1687 - eps(3)) && (off < 0.1687 + eps(3))) {  Tmin= 10; Tmax=19;     }
    else if( (tht > 0.2975 -eps(3)) && (tht < 0.2975 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 10; Tmax=19;     }
    else if( (tht > 0.4064 -eps(3)) && (tht < 0.4064 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 11; Tmax=19;     }
    else if( (tht > 0.4620 -eps(3)) && (tht < 0.4620 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 11; Tmax=19;     }
    else if( (tht > 0.3005 -eps(3)) && (tht < 0.3005 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 14; Tmax=22;     }
    else if( (tht > 0.4234 -eps(3)) && (tht < 0.4234 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 13; Tmax=19;     }
    else if( (tht > 0.4986 -eps(3)) && (tht < 0.4986 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 14; Tmax=20;     }
    else if( (tht > 0.5409 -eps(3)) && (tht < 0.5409 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 13; Tmax=20;     }
    else if( (tht > 0.2852 -eps(3)) && (tht < 0.2852 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 14; Tmax=22;     }
    else if( (tht > 0.4150 -eps(3)) && (tht < 0.4150 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 14; Tmax=22;     }
    else if( (tht > 0.5037 -eps(3)) && (tht < 0.5037 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 14; Tmax=20;     }
    else if( (tht > 0.5636 -eps(3)) && (tht < 0.5636 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 14; Tmax=20;     }
    else if( (tht > 0.5986 -eps(3)) && (tht < 0.5986 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 14; Tmax=20;     }
    else crash("In Get_virtual_ff_fit_interval I cannot find the kinematical point: off: "+to_string(off)+" , tht: "+to_string(tht));
  }


  else crash ("Unrecognized form factor type: "+W+" in Get_virtual_ff_fit_interval");




  return;
}


/* void Get_virtual_ff_fit_interval_v1(string W, double off, double tht, int &Tmin, int &Tmax) {


  if( W == "H1") {
    if( (tht > 0.2297 -eps(3)) && (tht < 0.2297 + eps(3)) && (off > 0.2003 - eps(3)) && (off < 0.2003 + eps(3))) { Tmin= 11; Tmax=15;     }
    else if( (tht > 0.2756 -eps(3)) && (tht < 0.2756 + eps(3)) && (off > 0.1687 - eps(3)) && (off < 0.1687 + eps(3))) {  Tmin= 13; Tmax=17;     }
    else if( (tht > 0.3603 -eps(3)) && (tht < 0.3603 + eps(3)) && (off > 0.1687 - eps(3)) && (off < 0.1687 + eps(3))) {  Tmin= 13; Tmax=17;     }
    else if( (tht > 0.2975 -eps(3)) && (tht < 0.2975 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 14; Tmax=18;     }
    else if( (tht > 0.4064 -eps(3)) && (tht < 0.4064 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 14; Tmax=17;     }
    else if( (tht > 0.4620 -eps(3)) && (tht < 0.4620 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 14; Tmax=17;     }
    else if( (tht > 0.3005 -eps(3)) && (tht < 0.3005 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 17; Tmax=21;     }
    else if( (tht > 0.4234 -eps(3)) && (tht < 0.4234 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 17; Tmax=20;     }
    else if( (tht > 0.4986 -eps(3)) && (tht < 0.4986 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 17; Tmax=20;     }
    else if( (tht > 0.5409 -eps(3)) && (tht < 0.5409 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 16; Tmax=20;     }
    else if( (tht > 0.2852 -eps(3)) && (tht < 0.2852 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 18; Tmax=21;     }
    else if( (tht > 0.4150 -eps(3)) && (tht < 0.4150 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 18; Tmax=21;     }
    else if( (tht > 0.5037 -eps(3)) && (tht < 0.5037 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 18; Tmax=21;     }
    else if( (tht > 0.5636 -eps(3)) && (tht < 0.5636 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 17; Tmax=20;     }
    else if( (tht > 0.5986 -eps(3)) && (tht < 0.5986 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 17; Tmax=20;     }
    else crash("In Get_virtual_ff_fit_interval I cannot find the kinematical point: off: "+to_string(off)+" , tht: "+to_string(tht));
  }

  else if (W == "H2") {
    if( (tht > 0.2297 -eps(3)) && (tht < 0.2297 + eps(3)) && (off > 0.2003 - eps(3)) && (off < 0.2003 + eps(3))) { Tmin= 11; Tmax=17;     }
    else if( (tht > 0.2756 -eps(3)) && (tht < 0.2756 + eps(3)) && (off > 0.1687 - eps(3)) && (off < 0.1687 + eps(3))) {  Tmin= 12; Tmax=20;     }
    else if( (tht > 0.3603 -eps(3)) && (tht < 0.3603 + eps(3)) && (off > 0.1687 - eps(3)) && (off < 0.1687 + eps(3))) {  Tmin= 14; Tmax=19;     }
    else if( (tht > 0.2975 -eps(3)) && (tht < 0.2975 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 15; Tmax=20;     }
    else if( (tht > 0.4064 -eps(3)) && (tht < 0.4064 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 15; Tmax=20;     }
    else if( (tht > 0.4620 -eps(3)) && (tht < 0.4620 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 15; Tmax=19;     }
    else if( (tht > 0.3005 -eps(3)) && (tht < 0.3005 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 16; Tmax=24;     }
    else if( (tht > 0.4234 -eps(3)) && (tht < 0.4234 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 15; Tmax=24;     }
    else if( (tht > 0.4986 -eps(3)) && (tht < 0.4986 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 15; Tmax=22;     }
    else if( (tht > 0.5409 -eps(3)) && (tht < 0.5409 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 15; Tmax=22;     }
    else if( (tht > 0.2852 -eps(3)) && (tht < 0.2852 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 16; Tmax=24;     }
    else if( (tht > 0.4150 -eps(3)) && (tht < 0.4150 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 16; Tmax=24;     }
    else if( (tht > 0.5037 -eps(3)) && (tht < 0.5037 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 16; Tmax=22;     }
    else if( (tht > 0.5636 -eps(3)) && (tht < 0.5636 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 16; Tmax=22;     }
    else if( (tht > 0.5986 -eps(3)) && (tht < 0.5986 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 15; Tmax=20;     }
    else crash("In Get_virtual_ff_fit_interval I cannot find the kinematical point: off: "+to_string(off)+" , tht: "+to_string(tht));


  }

  else if (W == "FA") {
    if( (tht > 0.2297 -eps(3)) && (tht < 0.2297 + eps(3)) && (off > 0.2003 - eps(3)) && (off < 0.2003 + eps(3))) { Tmin= 15; Tmax=18;     }
    else if( (tht > 0.2756 -eps(3)) && (tht < 0.2756 + eps(3)) && (off > 0.1687 - eps(3)) && (off < 0.1687 + eps(3))) {  Tmin= 15; Tmax=18;     }
    else if( (tht > 0.3603 -eps(3)) && (tht < 0.3603 + eps(3)) && (off > 0.1687 - eps(3)) && (off < 0.1687 + eps(3))) {  Tmin= 14; Tmax=17;     }//17 unsure
    else if( (tht > 0.2975 -eps(3)) && (tht < 0.2975 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 14; Tmax=18;     }
    else if( (tht > 0.4064 -eps(3)) && (tht < 0.4064 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 14; Tmax=17;     }//17 unsure
    else if( (tht > 0.4620 -eps(3)) && (tht < 0.4620 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 14; Tmax=16;     }//16 unsure
    else if( (tht > 0.3005 -eps(3)) && (tht < 0.3005 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 8; Tmax=12;     }
    else if( (tht > 0.4234 -eps(3)) && (tht < 0.4234 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 9; Tmax=13;     }
    else if( (tht > 0.4986 -eps(3)) && (tht < 0.4986 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 9; Tmax=13;     }
    else if( (tht > 0.5409 -eps(3)) && (tht < 0.5409 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 9; Tmax=13;     }
    else if( (tht > 0.2852 -eps(3)) && (tht < 0.2852 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 8; Tmax=13;     }
    else if( (tht > 0.4150 -eps(3)) && (tht < 0.4150 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 9; Tmax=13;     }
    else if( (tht > 0.5037 -eps(3)) && (tht < 0.5037 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 9; Tmax=13;     }
    else if( (tht > 0.5636 -eps(3)) && (tht < 0.5636 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 9; Tmax=13;     }
    else if( (tht > 0.5986 -eps(3)) && (tht < 0.5986 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 9; Tmax=13;     }
    else crash("In Get_virtual_ff_fit_interval I cannot find the kinematical point: off: "+to_string(off)+" , tht: "+to_string(tht));
   


  }

  else if (W =="FV") {
    if( (tht > 0.2297 -eps(3)) && (tht < 0.2297 + eps(3)) && (off > 0.2003 - eps(3)) && (off < 0.2003 + eps(3))) { Tmin= 11; Tmax=17;     }
    else if( (tht > 0.2756 -eps(3)) && (tht < 0.2756 + eps(3)) && (off > 0.1687 - eps(3)) && (off < 0.1687 + eps(3))) {  Tmin= 11; Tmax=19;     }
    else if( (tht > 0.3603 -eps(3)) && (tht < 0.3603 + eps(3)) && (off > 0.1687 - eps(3)) && (off < 0.1687 + eps(3))) {  Tmin= 13; Tmax=19;     }
    else if( (tht > 0.2975 -eps(3)) && (tht < 0.2975 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 13; Tmax=19;     }
    else if( (tht > 0.4064 -eps(3)) && (tht < 0.4064 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 13; Tmax=19;     }
    else if( (tht > 0.4620 -eps(3)) && (tht < 0.4620 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 13; Tmax=19;     }
    else if( (tht > 0.3005 -eps(3)) && (tht < 0.3005 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 15; Tmax=20;     }
    else if( (tht > 0.4234 -eps(3)) && (tht < 0.4234 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 15; Tmax=20;     }
    else if( (tht > 0.4986 -eps(3)) && (tht < 0.4986 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 15; Tmax=20;     }
    else if( (tht > 0.5409 -eps(3)) && (tht < 0.5409 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 15; Tmax=20;     }
    else if( (tht > 0.2852 -eps(3)) && (tht < 0.2852 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 15; Tmax=20;     }
    else if( (tht > 0.4150 -eps(3)) && (tht < 0.4150 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 15; Tmax=20;     }
    else if( (tht > 0.5037 -eps(3)) && (tht < 0.5037 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 15; Tmax=20;     }
    else if( (tht > 0.5636 -eps(3)) && (tht < 0.5636 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 15; Tmax=20;     }
    else if( (tht > 0.5986 -eps(3)) && (tht < 0.5986 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 15; Tmax=20;     }
    else crash("In Get_virtual_ff_fit_interval I cannot find the kinematical point: off: "+to_string(off)+" , tht: "+to_string(tht));
  }


  else crash ("Unrecognized form factor type: "+W+" in Get_virtual_ff_fit_interval");




  return;
}

*/

void Get_virtual_ff_fit_interval_v1(string W, double off, double tht, int &Tmin, int &Tmax) {


  if( W == "H1") {
    if( (tht > 0.2297 -eps(3)) && (tht < 0.2297 + eps(3)) && (off > 0.2003 - eps(3)) && (off < 0.2003 + eps(3))) { Tmin= 11; Tmax=15;     }
    else if( (tht > 0.2756 -eps(3)) && (tht < 0.2756 + eps(3)) && (off > 0.1687 - eps(3)) && (off < 0.1687 + eps(3))) {  Tmin= 13; Tmax=17;     }
    else if( (tht > 0.3603 -eps(3)) && (tht < 0.3603 + eps(3)) && (off > 0.1687 - eps(3)) && (off < 0.1687 + eps(3))) {  Tmin= 13; Tmax=17;     }
    else if( (tht > 0.2975 -eps(3)) && (tht < 0.2975 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 14; Tmax=17;     }
    else if( (tht > 0.4064 -eps(3)) && (tht < 0.4064 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 13; Tmax=17;     }
    else if( (tht > 0.4620 -eps(3)) && (tht < 0.4620 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 13; Tmax=17;     }
    else if( (tht > 0.3005 -eps(3)) && (tht < 0.3005 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 14; Tmax=20;     }
    else if( (tht > 0.4234 -eps(3)) && (tht < 0.4234 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 14; Tmax=20;     }
    else if( (tht > 0.4986 -eps(3)) && (tht < 0.4986 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 14; Tmax=20;     }
    else if( (tht > 0.5409 -eps(3)) && (tht < 0.5409 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 14; Tmax=20;     }
    else if( (tht > 0.2852 -eps(3)) && (tht < 0.2852 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 14; Tmax=20;     }
    else if( (tht > 0.4150 -eps(3)) && (tht < 0.4150 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 14; Tmax=21;     }
    else if( (tht > 0.5037 -eps(3)) && (tht < 0.5037 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 14; Tmax=20;     }
    else if( (tht > 0.5636 -eps(3)) && (tht < 0.5636 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 14; Tmax=20;     }
    else if( (tht > 0.5986 -eps(3)) && (tht < 0.5986 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 14; Tmax=20;     }
    else crash("In Get_virtual_ff_fit_interval I cannot find the kinematical point: off: "+to_string(off)+" , tht: "+to_string(tht));
  }

  else if (W == "H2") {
    if( (tht > 0.2297 -eps(3)) && (tht < 0.2297 + eps(3)) && (off > 0.2003 - eps(3)) && (off < 0.2003 + eps(3))) { Tmin= 11; Tmax=17;     }//horrible one
    else if( (tht > 0.2756 -eps(3)) && (tht < 0.2756 + eps(3)) && (off > 0.1687 - eps(3)) && (off < 0.1687 + eps(3))) {  Tmin= 12; Tmax=20;     }
    else if( (tht > 0.3603 -eps(3)) && (tht < 0.3603 + eps(3)) && (off > 0.1687 - eps(3)) && (off < 0.1687 + eps(3))) {  Tmin= 14; Tmax=18;     }
    else if( (tht > 0.2975 -eps(3)) && (tht < 0.2975 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 15; Tmax=20;     }
    else if( (tht > 0.4064 -eps(3)) && (tht < 0.4064 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 15; Tmax=19;     }
    else if( (tht > 0.4620 -eps(3)) && (tht < 0.4620 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 15; Tmax=19;     }
    else if( (tht > 0.3005 -eps(3)) && (tht < 0.3005 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 16; Tmax=24;     }//from here to the end all plateaux are "double" and we cover both
    else if( (tht > 0.4234 -eps(3)) && (tht < 0.4234 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 15; Tmax=24;     }
    else if( (tht > 0.4986 -eps(3)) && (tht < 0.4986 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 15; Tmax=22;     }
    else if( (tht > 0.5409 -eps(3)) && (tht < 0.5409 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 15; Tmax=22;     }
    else if( (tht > 0.2852 -eps(3)) && (tht < 0.2852 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 16; Tmax=24;     }
    else if( (tht > 0.4150 -eps(3)) && (tht < 0.4150 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 16; Tmax=22;     }
    else if( (tht > 0.5037 -eps(3)) && (tht < 0.5037 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 15; Tmax=22;     }
    else if( (tht > 0.5636 -eps(3)) && (tht < 0.5636 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 15; Tmax=22;     }
    else if( (tht > 0.5986 -eps(3)) && (tht < 0.5986 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 15; Tmax=22;     }
    else crash("In Get_virtual_ff_fit_interval I cannot find the kinematical point: off: "+to_string(off)+" , tht: "+to_string(tht));


  }

  else if (W == "FA") {
    if( (tht > 0.2297 -eps(3)) && (tht < 0.2297 + eps(3)) && (off > 0.2003 - eps(3)) && (off < 0.2003 + eps(3))) { Tmin= 15; Tmax=18;     }
    else if( (tht > 0.2756 -eps(3)) && (tht < 0.2756 + eps(3)) && (off > 0.1687 - eps(3)) && (off < 0.1687 + eps(3))) {  Tmin= 15; Tmax=18;     }
    else if( (tht > 0.3603 -eps(3)) && (tht < 0.3603 + eps(3)) && (off > 0.1687 - eps(3)) && (off < 0.1687 + eps(3))) {  Tmin= 14; Tmax=17;     }
    else if( (tht > 0.2975 -eps(3)) && (tht < 0.2975 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 15; Tmax=19;     }
    else if( (tht > 0.4064 -eps(3)) && (tht < 0.4064 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 14; Tmax=17;     }
    else if( (tht > 0.4620 -eps(3)) && (tht < 0.4620 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 13; Tmax=16;     }//short one
    else if( (tht > 0.3005 -eps(3)) && (tht < 0.3005 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 10; Tmax=13;     }//weird one
    else if( (tht > 0.4234 -eps(3)) && (tht < 0.4234 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 7; Tmax=13;     }//long one ?
    else if( (tht > 0.4986 -eps(3)) && (tht < 0.4986 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 8; Tmax=13;     }
    else if( (tht > 0.5409 -eps(3)) && (tht < 0.5409 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 8; Tmax=13;     }
    else if( (tht > 0.2852 -eps(3)) && (tht < 0.2852 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 10; Tmax=13;     }//weird one
    else if( (tht > 0.4150 -eps(3)) && (tht < 0.4150 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 10; Tmax=14;     }
    else if( (tht > 0.5037 -eps(3)) && (tht < 0.5037 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 8; Tmax=13;     }
    else if( (tht > 0.5636 -eps(3)) && (tht < 0.5636 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 8; Tmax=13;     }
    else if( (tht > 0.5986 -eps(3)) && (tht < 0.5986 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 8; Tmax=13;     }
    else crash("In Get_virtual_ff_fit_interval I cannot find the kinematical point: off: "+to_string(off)+" , tht: "+to_string(tht));
   


  }

  else if (W =="FV") {
    if( (tht > 0.2297 -eps(3)) && (tht < 0.2297 + eps(3)) && (off > 0.2003 - eps(3)) && (off < 0.2003 + eps(3))) { Tmin= 11; Tmax=17;     }
    else if( (tht > 0.2756 -eps(3)) && (tht < 0.2756 + eps(3)) && (off > 0.1687 - eps(3)) && (off < 0.1687 + eps(3))) {  Tmin= 11; Tmax=19;     }
    else if( (tht > 0.3603 -eps(3)) && (tht < 0.3603 + eps(3)) && (off > 0.1687 - eps(3)) && (off < 0.1687 + eps(3))) {  Tmin= 13; Tmax=19;     }
    else if( (tht > 0.2975 -eps(3)) && (tht < 0.2975 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 13; Tmax=19;     }
    else if( (tht > 0.4064 -eps(3)) && (tht < 0.4064 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 13; Tmax=19;     }
    else if( (tht > 0.4620 -eps(3)) && (tht < 0.4620 + eps(3)) && (off > 0.1370 - eps(3)) && (off < 0.1370 + eps(3))) {  Tmin= 13; Tmax=19;     }
    else if( (tht > 0.3005 -eps(3)) && (tht < 0.3005 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 15; Tmax=20;     }
    else if( (tht > 0.4234 -eps(3)) && (tht < 0.4234 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 15; Tmax=20;     }
    else if( (tht > 0.4986 -eps(3)) && (tht < 0.4986 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 15; Tmax=20;     }
    else if( (tht > 0.5409 -eps(3)) && (tht < 0.5409 + eps(3)) && (off > 0.1054 - eps(3)) && (off < 0.1054 + eps(3))) {  Tmin= 15; Tmax=20;     }
    else if( (tht > 0.2852 -eps(3)) && (tht < 0.2852 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 15; Tmax=20;     }
    else if( (tht > 0.4150 -eps(3)) && (tht < 0.4150 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 15; Tmax=20;     }
    else if( (tht > 0.5037 -eps(3)) && (tht < 0.5037 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 15; Tmax=20;     }
    else if( (tht > 0.5636 -eps(3)) && (tht < 0.5636 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 15; Tmax=20;     }
    else if( (tht > 0.5986 -eps(3)) && (tht < 0.5986 + eps(3)) && (off > 0.0738 - eps(3)) && (off < 0.0738 + eps(3))) {  Tmin= 15; Tmax=20;     }
    else crash("In Get_virtual_ff_fit_interval I cannot find the kinematical point: off: "+to_string(off)+" , tht: "+to_string(tht));
  }


  else crash ("Unrecognized form factor type: "+W+" in Get_virtual_ff_fit_interval");




  return;
}
