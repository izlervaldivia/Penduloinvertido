function xdot=edsonj(u1,u2,u3,u4,u5)
// This is nonlinear EDSON-J model

// Load the parameters
exec('edsonjParameters.sce', -1);
  m=0.017;
  M=0.696;
  l=0.3;
  g=9.81;
// state variables

//u1;		//angulo
//u2;     //velocidad angular
//u3;      //posicion
//u4;      // velocidad


// control variables
//u5	

//
xdot=[  u2;
       (+l*m*cos(u1)*sin(u1)*u2^2 + u*cos(u1)-g*sin(u1)*(M+m))/((m*l*cos(u1)^2)-((M+m)*l));
       u4;
       (m*l*u2^2*sin(u1)-g*m*cos(u1)*sin(u1)+u5)/(M+m+(m*cos(u1)^2))];

endfunction
