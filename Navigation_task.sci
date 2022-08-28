clear;
clc;

//1k = 0.1s 1800k=3min
//0.01 = 1cm 1=1m

function [dxdt]=DynamicModel(t,c)
    // c(1,:) x position
    // c(2,:) y position
    // c(3,:) movement
    // c(4,:) orientation
    // c(5,:) speed
    // c(6,:) angular speed
    
    // parameters
    m = .38;      //mass
    Ip = 0.005; //inertia moment
    d = 0.02;      //distance centroide to axis wheel 
    r = 0.03;      //radio wheel 
    R = 0.05;      //distance wheel-center 
    
    
    M = [m 0 ; 0 (Ip+m*d^2)];         //inertial matrix
    C = [-m*d*c(6)^2 ; m*d*c(5)*c(6)] //matrix of Coriolis and centrifugal forces 
    B = [1/r 1/r ; R/r -R/r];         //matriz conversion  torque-wheel-movil force
    A = [r/2 r/2 ; r/(2*R) -r/(2*R)];
    Ts = [.434 0 ; 0 .434];
    Ks = [2.745 0 ; 0 2.745];
    Kl = [1460.2705 0 ; 0 1460.2705];
       
    dxdt=[[cos(c(4)) -d*sin(c(4)); sin(c(4)) d*cos(c(4))]*[c(5); c(6)];
    c(5);
    c(6);
    inv(M + B*inv(Kl)*Ts*inv(A)) * (B*inv(Kl)*Ks*u - (C + B*inv(Kl)*inv(A)*[c(5);c(6)]))];
endfunction

function [c,t]= Mov(u,ci) //movil dynamic
    //initial variables
    //global tau; //torque in every wheel a vector default tau=[1;1];
    //tau=Torque;
    ti=0;//initial time
    tf=1;//final time
    tspan=[ti:0.1:tf];//tme vector
    
//    errcatch(98,"pause");
    c=ode(ci,ti,tspan,DynamicModel);
    t=tspan;
endfunction


function [res]=SimMov2(Tv,dR,dO,dA)
    Tm= 180;
    Report=zeros(Tm,Tv,3);//reporte
    AlphaR=zeros(Tm,Tv,1);//reporte alpha
    res=zeros(3,5); //respultados centro de masa y dispersion
    C=zeros(Tv,6); //estados del vehiculo
    rm = rand(1,1,'normal') * 0.01;
    uO =2; uR =1; uA =3.2; uI =3.2 ;//velocidades
    ud = zeros(Tv,2)+uO;//velocidad deseada
    limsx= 24; //limites
    limsy= 4;
    //dI = (lims*0.5); // distancia de influencia
    dI = (limsx*0.8);
    N1=(Tm/2);
    
    figure()
    plot(limsx,limsy,'k');
    f=get("current_figure");
    f.background=8;
    xlabel('x-axis')
    ylabel('y-axis')

   //Inicializar variables
    for v=1:Tv
        C(v,1)=ceil(limsx*0.2*rand());       //x position
        C(v,2)=ceil(limsy*rand());       //y position
        C(v,3)=0;            //movement made
        C(v,4)=2*%pi*rand(); //orientation
        C(v,5)=0;            //Speed
        C(v,6)=0;            //angular speed 
    end // end v111
    
    for m=1:Tm
        
        //posicion de influencia
        if m<Tm then
            inf =[limsx*0.9,limsy*0.5];
        end
        // valor de influencia activa y pasiva
        if m<(Tm*0.05) then
            dInf = 0;
        else
            dInf = dI;
        end
        
        for v=1:Tv   
            [cs,t] = Mov(ud(v,:)',C(v,:)');
            [r, c]=size(cs);
            C(v,:)=cs(:,c)';//next initial condition, actual condition
            ud(v,:) = zeros(1,2) + uO;
            
            //this avoids an infinite increment of radians
            if C(v,4)>(2*%pi) then
                C(v,4)=C(v,4)-(2*%pi);
            end
            if C(v,4)<0 then
                C(v,4)=C(v,4)+(2*%pi);
            end
            
            //Limits
            turnLIM = (%pi/2);
            if C(v,1)>(limsx-(limsx*0.05)) then
                if (0<((C(v,4)*180)/%pi) & ((C(v,4)*180)/%pi)<90) then
                    C(v,4)=C(v,4)+(0.5*rand()+turnLIM);
                end
                
                if (270<((C(v,4)*180)/%pi) & ((C(v,4)*180)/%pi)<360) then
                    C(v,4)=C(v,4)+(0.5*rand()-turnLIM);
                end
            end
            if C(v,2)>(limsy)-(limsy*0.05) then
                if (90<((C(v,4)*180)/%pi) & ((C(v,4)*180)/%pi)<180) then
                    C(v,4)=C(v,4)+(0.5*rand()+turnLIM);
                end
                if (0<((C(v,4)*180)/%pi) & ((C(v,4)*180)/%pi)<90) then
                    C(v,4)=C(v,4)+(0.5*rand()-turnLIM);
                end
            end
            if C(v,1)<(limsx*0.05) then
                if (180<((C(v,4)*180)/%pi) & ((C(v,4)*180)/%pi)<270) then
                    C(v,4)=C(v,4)+(0.5*rand()+turnLIM);
                end
                if (90<((C(v,4)*180)/%pi) & ((C(v,4)*180)/%pi)<180) then
                    C(v,4)=C(v,4)+(0.5*rand()-turnLIM);
                end
            end 
            if C(v,2)<(limsy*0.05) then
                if (270<((C(v,4)*180)/%pi) & ((C(v,4)*180)/%pi)<360) then
                    C(v,4)=C(v,4)+(0.5*rand()+turnLIM);
                end
                if (180<((C(v,4)*180)/%pi) & ((C(v,4)*180)/%pi)<270) then
                    C(v,4)=C(v,4)+(0.5*rand()-turnLIM);
                end
            end
             
            //Elements detected
            DetR=0; DetA=0; DetI=0; DetO=0;
            ElemR = []; ElemA = [];
            
            //Rango de detección para repulsión, orientación, atracción e influencia.
            rangR=1.74533; rangA= 0.5235988; rangI=2.79253; 

            for z=1:Tv
                // Angulo del individuo con respecto a otros miembros del enjambre
               angD = atan(C(z,2)-C(v,2),C(z,1)-C(v,1))
               if angD<0 then
                   angD=angD+(2*%pi);
               end
               if angD > (2*%pi) then
                   angD=angD-(2*%pi);
               end
               
               //Calculo de angulos rep y atrac
               Beta = angD - C(v,4);
               if  Beta < 0 then
                   Beta = Beta + (2*%pi);
               end
               
               Gama = C(v,4) - angD;
               if  Gama < 0 then
                   Gama = Gama + (2*%pi);
               end
               
               if Gama<Beta then
                   Delta = Gama;
               else
                   Delta = Beta;
               end
               
               //Calculo distancia repulsion
               if Delta<rangR then
                   reps=1
                   if v<>z then //it must not be the same
                       dr=sqrt((C(v,1)-C(z,1))^2+(C(v,2)-C(z,2))^2);
                   else
                       dr=1000000;
                   end
               else
                   dr=1000000;
               end
               
               //Calculo distancia atraccion
               if Delta<rangA then
                   if v<>z then //it must not be the same
                       da=sqrt((C(v,1)-C(z,1))^2+(C(v,2)-C(z,2))^2);
                   else
                       da=1000000;
                   end
               else
                   da=1000000;
               end
                
           //Numero de individuos detectados en el radio de atraccion y repulsion
                if  dr<dR then
                    ElemR = [ElemR ; angD];
                    DetR=DetR+1; //A vehicle is detected
                end
                
                if  da>dO & da<dA & DetR==0 then
                    ElemA = [ElemA ; angD];
                    DetA=DetA+1; //A vehicle is detected
                end
                
                d=sqrt((C(v,1)-C(z,1))^2+(C(v,2)-C(z,2))^2);
                if d>dR & d<dO & DetR==0 & DetA ==0 then
                    DetO=DetO+1; //A vehicle is detected
                end
            end //of z

            //Ángulo de influencia
            angI = atan(inf(2)-C(v,2),inf(1)-C(v,1));
            if angI<0 then
               angI=angI+(2*%pi);
            end
            if angI > (2*%pi) then
               angI=angI-(2*%pi);
            end
            
               //Calculo de angulos rep y atrac
               BetaI = angI - C(v,4);
               if  BetaI < 0 then
                   BetaI = BetaI + (2*%pi);
               end
               
               GamaI = C(v,4) - angI;
               if  GamaI < 0 then
                   GamaI = GamaI + (2*%pi);
               end
               
               if GamaI<BetaI then
                   DeltaI = GamaI;
               else
                   DeltaI = BetaI;
               end

               //ruido para valor de distancia en influencia (distribución estandar) del 5%
               suma = 0;
               resta = 0;
               for i=1 : 12
                   suma =  suma + rand(1,1,'normal');
               end
               for j=1 : 6
                   resta = resta + rand(1,1,'normal');
               end
               z= (suma -resta) * .05;
            
           //Distancia de influencia calculada            
            if Delta<rangI then
                di=sqrt((inf(1)-C(v,1))^2+(inf(2)-C(v,2))^2) + z;
            else
                di = 100000;
            end
            
            //Promedio de elementos detectados
            PromElemR = mean(ElemR)
            PromElemA = mean(ElemA)
            
            //REPULSION
            if DetR > 0 then
               ud(v,:) = zeros(1,2) + uR + rm;
               if PromElemR < %pi then
                   C(v,4)=PromElemR + %pi;
               else
                   C(v,4)=PromElemR - %pi;
               end
            end
           
            //INFLUENCIA
            if di<dInf & DetR==0 then
                ud(v,:) = zeros(1,2) + uI + rm;
                C(v,4) = (angI*0.6) + (C(v,4)*0.4)
                //C(v,4) = angI
                DetI==1;
            end
           
            //ATRACCIÓN
            if DetA>0 & DetR==0 & DetI==0 then
                ud(v,:) = zeros(1,2) + uA + rm;
                C(v,4)=PromElemA;
            end
            
            //ORIENTACIÓN
            if DetO>0 & DetR==0 & DetA==0 & DetI==0 then
                ud(v,:) = zeros(1,2) + uO + rm;
                C(v,4) = C(v,4);
            end
           
            //fuera de rango
            if DetO==0 & DetR==0 & DetA==0 & DetI==0 then
                ud(v,:) = zeros(1,2) + uO + rm;
                C(v,4) = C(v,4);
            end

           //calculo grados para orientacion de flechas o robots
           phi = C(v,4);
           angles=(phi*180)/%pi;
           alfa = C(v,4) - (0.5*%pi);
           
            if alfa>(2*%pi) then
                alfa=alfa-(2*%pi);
            end
            if alfa<0 then
                alfa=alfa+(2*%pi);
            end
           
           //Reports
           Report(m,v,1) = C(v,1);
           Report(m,v,2) = C(v,2);
           Report(m,v,3) = angles;
           Report(m,v,4) = C(v,5);
          
          //reporte orientacion de flechas
           AlphaR(m,v,1) = alfa

         end //of v 
        //Plot vehicles positions
        a = gca();delete(a.children);
        for v=1:Tv 
            xf(v)=C(v,1)-0.01*sin(AlphaR(m,v,1));
            yf(v)=C(v,2)+0.01*cos(AlphaR(m,v,1));
        end
         nx=[C(:,1)';xf'];
         ny=[C(:,2)';yf'];
         xarrows(nx,ny,(Tv/1)) 
        set(gca(),"auto_clear","off")
        h_compound = gce();
        if m>(Tm*0.05) then
            plot(inf(1),inf(2),'y.');
            h_compound = gce();
            h_compound.children.mark_size = 10;
        end

        if m==N1 then
            pnx1=nx;
            pny1=ny;
        end
        
    end //of m
     
     figure()
     f=get("current_figure");
     f.background=8;
     plot(limsx,limsy);
     plot(0,0);
        xlabel('x-axis')
        ylabel('y-axis')
     set(gca(),"auto_clear","off")
     xarrows(pnx1,pny1,(Tv/1));
     cx = sum(Report(N1,:,1))/Tv;
     cy = sum(Report(N1,:,2))/Tv;
     dy=sqrt(sum((Report(N1,:,2)-cy).^2)/Tv);
     dx=sqrt(sum((Report(N1,:,1)-cx).^2)/Tv);
     area=%pi*dx*dy
     res(1,1)=cx;
     res(1,2)=cy;
     res(1,3)=dx;
     res(1,4)=dy;
     res(1,5)=area;
     set(gca(),"auto_clear","off")
     plot(cx,cy,'rx')
     xset("color",0)
     xarc(cx-(dx*1.3),cy+(dy*1.3),(dx*1.3)*2,(dy*1.3)*2,0,360*64)
     set(gca(),"auto_clear","off")
     plot(inf(1),inf(2),'y.');
     h_compound = gce();
     h_compound.children.mark_size = 10;
     set(gca(),"auto_clear","on")
endfunction

function [promres,ress]=simulacion(Tv,dR,dO,dA)
    for i=1:3
        [res]=SimMov2(Tv,dR,dO,dA)
        for j=1:1
            for k=1:5
                ress(i,j,k)=res(j,k);
            end
        end
    end
    for i=1:1
        for j=1:5
           promres(i,j)= sum(ress(:,i,j))/3;
        end
    end


endfunction
