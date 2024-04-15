%LAMMPS金刚石纳米柱建模文件 (竖直方向） 
%拉伸方向 011

close all
clear all

name='diamond_data.in';
l=1.54; %C-C bond length
RG=10;  %cutting radius
H=30;   %cutting height
%% 

diab=[0,0,0];
latticeX=20;
latticeY=20;   
latticeZ=300; 


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
dia=zeros(p,4);
p=1;



%% establish 3-D diamond structure 
for k=1:latticeZ    
    for j=1:latticeY
        for i=1:latticeX
            dia(p,1)=diab(1)+(i-1)*dx;
            dia(p,2)=diab(2)+(j-1)*dy;
            dia(p,3)=diab(3)+(k-1)*dz;
            if mod(i,2)==1
               if mod(k,2)==1
               else
                   dia(p,2)=dia(p,2)-c;
               end
            else 
               if mod(k,2)==1
                   dia(p,2)=dia(p,2)+c;
               else 
                   dia(p,2)=dia(p,2)-2*c;
               end
            end
            p=p+1;
        end
    end
end

p=p-1;  %total nomimal atoms number

Xdia=[min(dia(:,1)),max(dia(:,1))]; % boundary of the diamond
Ydia=[min(dia(:,2)),max(dia(:,2))];
Zdia=[min(dia(:,3)),max(dia(:,3))];

RW=Xdia(2)-Xdia(1); % width and length
RL=Ydia(2)-Ydia(1);

XGmid=mean(dia(:,1));
YGmid=mean(dia(:,2));
ZGmid=mean(dia(:,3));

ng=p;

ii=0;                   %切半径
diaR=zeros(p,4);  
for i=1:p
    xdis=abs(dia(i,1)-XGmid);
    ydis=abs(dia(i,2)-YGmid);
    dist=sqrt(xdis^2+ydis^2);
    if dist<=RG
        ii=ii+1;
        diaR(ii,:)=dia(i,:);
    end
    
end
ng=ii;
diaR(ii+1:p,:)=[];
dia=diaR;
dia(:,4)=1;
ii=0;         
mm=max(dia(:,2))
mini=min(dia(:,2))           %切标记
for i=1:ng
 
    if dia(i,2)==mm && dia(i,1)>XGmid
        dia(i,4)=4;
    end
    if dia(i,2)==mini && dia(i,1)>XGmid
        dia(i,4)=5;
    end
end


diaH=zeros(p,4);
H=H/2;               %切上下端
ii=0;
for i=1:ng
    h=abs(dia(i,3)-ZGmid);
    if h<H
        ii=ii+1;
        diaH(ii,:)=dia(i,:);
    end
    
end
ng=ii;
diaH(ii+1:p,:)=[];
dia=diaH;
            
%define the top boundary

BOT=dia(1,3)+dz+0.1;
UP=dia(ii,3)-dz-0.1;
for i=1:ng
    h=dia(i,3);
    if h>UP
        
        dia(i,4)=2;
    end
    
end

for i=1:ng
    h=dia(i,3);
    if h<BOT
        
        dia(i,4)=3;
    end
    
end



Xdia=[min(dia(:,1)),max(dia(:,1))]; % boundary of the diamond
Ydia=[min(dia(:,2)),max(dia(:,2))];
Zdia=[min(dia(:,3)),max(dia(:,3))];

hh=Zdia(2)-Zdia(1);
kk=floor(hh/dz);
nn=mod(kk,2);
if nn==1
    Zdia(1)=Zdia(1)-dz;
else
    Zdia(2)=Zdia(2)+dz;
end
    

%% OUTPUT
fid=fopen(name,'wt'); %Creat a file to write
fprintf(fid,'#Lammps data file\n\n');
fprintf(fid,'%5d atoms\n',ng); % totoal atoms number
fprintf(fid,'%5d atom types\n\n',5); % types number

% bondaries
fprintf(fid,'%8.4f %8.4f xlo xhi\n',Xdia(1)-10,Xdia(2)+10);
fprintf(fid,'%8.4f %8.4f ylo yhi\n',Ydia(1)-10,Ydia(2)+10);
fprintf(fid,'%8.4f %8.4f zlo zhi\n\n',Zdia(1),Zdia(2));

fprintf(fid,'Masses\n\n');
fprintf(fid,' 1 12\n'); 
fprintf(fid,' 2 12\n'); 
fprintf(fid,' 3 12\n'); 
fprintf(fid,' 4 12\n'); 
fprintf(fid,' 5 12\n\n'); 

fprintf(fid,'Atoms\n\n');

for i=1:ng
    out=i
    fprintf(fid,'%5d  %1d  %7.3f  %7.3f  %7.3f\n',i,dia(i,4),dia(i,1),dia(i,2),dia(i,3));
 
    
end

fclose(fid);
ok='OK'
            

            
         
            
        
        
        
        
        
        
        
        
        