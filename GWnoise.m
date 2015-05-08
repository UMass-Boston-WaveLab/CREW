function [ Hn ] = GWnoise( input_args )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
   XX=0
   for i = 1 to N
      U = rand()
      XX = XX + U
   end       
do		            % if U1==0, then the log() below will blow up
          U1=rand()   
while U1==0
   U2=rand()
	
   XX=sqrt(-2 * log(U1)) * cos(2pi * U2)
   YY=sqrt(-2 * log(U1)) * sin(2pi * U2)
  
do
      U1=rand()            /* U1=[0,1] */
      U2=rand()            /* U2=[0,1] */
      V1=2 * U1 -1            /* V1=[-1,1] */
      V2=2 * U2 - 1           /* V2=[-1,1] */
      S=V1 * V1 + V2 * V2
   while S >=1
	
   XX=sqrt(-2 * log(S) / S) * V1
   YY=sqrt(-2 * log(S) / S) * V2
end

