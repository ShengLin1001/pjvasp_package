
RG=r;  %cutting radius
H=h;   %height of the cylinder
l=1.54; %C-C bond length 
Grab=[0,0,0];
latticeX=200;
latticeY=200;   
latticeZ=400; %assume height extends in z direction
%% 

aa = 4.08196932030572;
bb = sqrt(3.0);
cc = 1.666;
hcp_a = aa * sqrt(2)/2;
base = hcp * [1.0 0 0;0 0 bb;0 0 cc]

%% 1


lattice custom ${hcp_a} a1 1.0 0 0 a2 0.0 ${bb} 0 a3 0 0 ${cc} &
                        basis 0.1   0.1             0.1                      &
                        basis 0.6   0.6           0.1                      &
                        basis 0.6   0.93333333333 0.6                    &
                        basis 0.1   0.43333333333 0.6 

p=1; %atom account
thita=(109+28/60)/180*pi;
a=cos(thita/2);  %cos(109.28/2)
b=sin(thita/2);  %sin(109.29/2)
c=a*l;
d=b*l;
dx=d;
dy=4*c;
dz=d;

p=latticeX*latticeY*latticeZ;  %total number 
Gra=zeros(p,4);
p=1;



%establish 3-D diamond structure 
for k=1:latticeZ    
    for j=1:latticeY
        for i=1:latticeX
            Gra(p,1)=Grab(1)+(i-1)*dx;
            Gra(p,2)=Grab(2)+(j-1)*dy;
            Gra(p,3)=Grab(3)+(k-1)*dz;
            if mod(i,2)==1
               if mod(k,2)==1
               else
                   Gra(p,2)=Gra(p,2)-c;
               end
            else 
               if mod(k,2)==1
                   Gra(p,2)=Gra(p,2)+c;
               else 
                   Gra(p,2)=Gra(p,2)-2*c;
               end
            end
            p=p+1;
        end
    end
end

p=p-1;  %total nomimal atoms number

XGra=[min(Gra(:,1)),max(Gra(:,1))]; % boundary of the diamond
YGra=[min(Gra(:,2)),max(Gra(:,2))];
ZGra=[min(Gra(:,3)),max(Gra(:,3))];

RW=XGra(2)-XGra(1); % width and length
RL=YGra(2)-YGra(1);

XGmid=mean(Gra(:,1));
YGmid=mean(Gra(:,2));
ZGmid=mean(Gra(:,3));

ng=p;

ii=0;                   %cutting the cylinder
GraR=zeros(p,4);  
for i=1:p
    xdis=abs(Gra(i,1)-XGmid);
    ydis=abs(Gra(i,2)-YGmid);
    dist=sqrt(xdis^2+ydis^2);
    if dist<=RG
        ii=ii+1;
        GraR(ii,:)=Gra(i,:);
    end
    
end
ng=ii;
GraR(ii+1:p,:)=[];
Gra=GraR;

GraH=zeros(p,4);   
H=H/2;               %cut by cylinder 
ii=0;
for i=1:ng
    h=abs(Gra(i,3)-ZGmid);
    if h<H
        ii=ii+1;
        GraH(ii,:)=Gra(i,:);
    end
    
end
ng=ii;
GraH(ii+1:p,:)=[];
GraH(:,4)=1;
Gra=GraH;
            
%define the top boundary

BOT=Gra(1,3)+dz+0.1;
UP=Gra(ii,3)-dz-0.1;
for i=1:ng
    h=Gra(i,3);
    if h>UP
        
        Gra(i,4)=2;
    end
    
end

for i=1:ng
    h=Gra(i,3);
    if h<BOT
        
        Gra(i,4)=3;
    end
    
end



XGra=[min(Gra(:,1)),max(Gra(:,1))]; % boundary of the diamond
YGra=[min(Gra(:,2)),max(Gra(:,2))];
ZGra=[min(Gra(:,3)),max(Gra(:,3))];

hh=ZGra(2)-ZGra(1);
kk=floor(hh/dz);
nn=mod(kk,2);
if nn==1
    ZGra(1)=ZGra(1)-dz;
else
    ZGra(2)=ZGra(2)+dz;
end
    


%OUTPUT
fid=fopen('diamond_datarh.in','wt'); %Creat a file to write
fprintf(fid,'#Lammps data file\n\n');
fprintf(fid,'%5d atoms\n',ng); % totoal atoms number
fprintf(fid,'%5d atom types\n\n',3); % types number

% bondaries
fprintf(fid,'%8.4f %8.4f xlo xhi\n',XGra(1)-10,XGra(2)+10);
fprintf(fid,'%8.4f %8.4f ylo yhi\n',YGra(1)-10,YGra(2)+10);
fprintf(fid,'%8.4f %8.4f zlo zhi\n\n',ZGra(1),ZGra(2));

fprintf(fid,'Masses\n\n');
fprintf(fid,' 1 12\n'); 
fprintf(fid,' 2 12\n'); 
fprintf(fid,' 3 12\n\n'); 

fprintf(fid,'Atoms\n\n');

for i=1:ng
    out=i
    fprintf(fid,'%5d  %1d  %7.3f  %7.3f  %7.3f\n',i,Gra(i,4),Gra(i,1),Gra(i,2),Gra(i,3));
 
    
end

fclose(fid);
ok='OK'
end