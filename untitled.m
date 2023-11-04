function varargout = untitled(varargin)
% UNTITLED MATLAB code for untitled.fig
%      UNTITLED, by itself, creates a new UNTITLED or raises the existing
%      singleton*.
%
%      H = UNTITLED returns the handle to a new UNTITLED or the handle to
%      the existing singleton*.
%
%      UNTITLED('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UNTITLED.M with the given input arguments.
%
%      UNTITLED('Property','Value',...) creates a new UNTITLED or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before untitled_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to untitled_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help untitled

% Last Modified by GUIDE v2.5 16-Feb-2023 19:29:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @untitled_OpeningFcn, ...
                   'gui_OutputFcn',  @untitled_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before untitled is made visible.
function untitled_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to untitled (see VARARGIN)

% Choose default command line output for untitled
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes untitled wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = untitled_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function fichier_Callback(hObject, eventdata, handles)
% hObject    handle to fichier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Bruit_Callback(hObject, eventdata, handles)
% hObject    handle to Bruit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------

% --------------------------------------------------------------------
function Filtre_Frequentiel_Callback(hObject, eventdata, handles)
% hObject    handle to Filtre_Frequentiel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function FPBI_Callback(hObject, eventdata, handles)
% hObject    handle to FPBI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function FPBB_Callback(hObject, eventdata, handles)
% hObject    handle to FPBB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


I = handles.courant_data;
%I = imread('eight.tif');

F=fftshift(fft2(I)); 
%calcul de la taille de l'image; 
M=size(F,1); 
N=size(F,2); 
P=size(F,3);

H0=zeros(M,N); 
D0=1; 
M2=round(M/2); 
N2=round(N/2); 
H0(M2-D0:M2+D0,N2-D0:N2+D0)=1; 

n=3; 

for i=1:M 
for j=1:N 
H(i,j)=1/(1+(H0(i,j)/D0)^(2*n)); 
G(i,j)=F(i,j)*H(i,j); 
end 
end 

g=ifft2(G); 

%subplot(1,2,1);imshow(I);title('image originale'); 
%subplot(1,2,2);
imshow(abs(g),[0,255]);%title('image filtrée'); 

% --------------------------------------------------------------------
function FPHB_Callback(hObject, eventdata, handles)
% hObject    handle to FPHB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I=handles.courant_data;

F=fftshift(fft2(I)); 
%calcul de la taille de l'image; 
M=size(F,1); 
N=size(F,2); 
P=size(F,3); 

H1=ones(M,N); 
D0=1; 
M2=round(M/2); 
N2=round(N/2); 
H1(M2-D0:M2+D0,N2-D0:N2+D0)=0; 

n=3; 

for i=1:M 
for j=1:N 
H(i,j)=1/(1+(H1(i,j)/D0)^(2*n)); 
G(i,j)=F(i,j)*H(i,j); 
end 
end 

g=ifft2(G); 

%subplot(1,2,1);imshow(I);title('image originale'); 
%subplot(1,2,2);
imshow(255-abs(g),[0,255]);%title('image filtrée');


% --------------------------------------------------------------------
function Gaussien_Callback(hObject, eventdata, handles)
% hObject    handle to Gaussien (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imageO=handles.courant_data;
Imgb =imnoise(imageO,'gaussian', 0,0.01);
%Imgb=handles.Imgb;
axes(handles.imgT)
imshow(Imgb);
handles.courant_data=Imgb;
%subImage(handles.Imgb);
handles.output=hObject;
guidata(hObject,handles);
title('image Bruite Gauss');


% --------------------------------------------------------------------
function puovre_et_sel_Callback(hObject, eventdata, handles)
% hObject    handle to puovre_et_sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imageO=handles.courant_data;
Imgb =imnoise(imageO,'Salt & Pepper', 0.01);
%Imgb=handles.Imgb;
axes(handles.imgT)
%
imshow(Imgb);
handles.courant_data=Imgb;
%subImage(handles.Imgb);
handles.output=hObject;
guidata(hObject,handles);
title('image Bruite Pov et Sel');



% --------------------------------------------------------------------

function Ouvrir_Callback(hObject, eventdata, handles)
% hObject    handle to Ouvrir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file,path] = uigetfile('*.*');
handles.ima = imread(sprintf('%s',path,file));
axes(handles.imgO)
handles.courant_data = handles.ima;
subimage(handles.courant_data);
title('image originale');
axes(handles.imgT)
% handles.ima_traite = 0;
% subimage(handles.ima_traite);

subimage(handles.courant_data);

handles.output = hObject;
guidata(hObject, handles);


% --------------------------------------------------------------------
function Enregistrer_Callback(hObject, eventdata, handles)
% hObject    handle to Enregistrer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = handles.courant_data;
[file,path] = uiputfile('*.png','Enregistrer Votre Image ...');
imwrite(image, sprintf('%s',path,file),'png');

% --------------------------------------------------------------------
function Quitter_Callback(hObject, eventdata, handles)
% hObject    handle to Quitter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1)

% --------------------------------------------------------------------
function Median_Callback(hObject, eventdata, handles)
% hObject    handle to Median (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
image=double(image);
[n,m]=size(image);
img=image;
for i=2:n-1
    for j=2:m-1
       fenetre=image(i-1:i+1,j-1:j+1);
       v=[fenetre(1,:) fenetre(2,:) fenetre(3,:)];
       sort(v);
       a=median(v);
       img(i,j)=a;
    end
end
b=uint8(img);
handles.ima_traite = b;
axes(handles.imgT);
subimage(b);
handles.output = hObject;
guidata(hObject, handles);
title('image filtre Median');
% --------------------------------------------------------------------


% --------------------------------------------------------------------
function Invesrtion_des_couleurs_Callback(hObject, eventdata, handles)
% hObject    handle to Invesrtion_des_couleurs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
%[n,m]=size(image);
image = double(image);
[l c]=size(image);
image = double(image);
v=image;
for i=1:l
   for j=1:c
     v(i,j)=-double(image(i,j))+255;
    end
 end 

v=uint8(v); 
axes(handles.imgT);
subimage(v);
title('image_inverse');


% --------------------------------------------------------------------
function Luminosite_Callback(hObject, eventdata, handles)
% hObject    handle to Luminosite (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[l c]=size(image);
image = double(image);
v=image;
for i=1:l
    for j=1:c
        pix=image(i,j)+50;
         if(pix>255)
            pix=255;
         else if (pix<0)
                pix=0;
             
              end 
          end
       v(i,j)=pix;    
    end
end  
v=uint8(v); 
axes(handles.imgT);
subimage(v);
title('image_luminosite');


% --------------------------------------------------------------------
function Niveau_de_gris_Callback(hObject, eventdata, handles)
% hObject    handle to Niveau_de_gris (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ima=handles.courant_data;
d = length(size(ima));
if d==3
    imagray=rgb2gray(ima);
elseif d==2
   imagray=ima;
end
axes(handles.imgT);
subimage(imagray);
title('image_Niveau_de_gris');
% --------------------------------------------------------------------
function menuHistog_Callback(hObject, eventdata, handles)
% hObject    handle to menuHistog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img = handles.courant_data;
%I=rgb2gray(img);
        d = length(size(img));
        if d==3
            I = rgb2gray(img);
        elseif d==2
            I = img
        end
axes(handles.imgO);
subimage(I);


[nl nc]=size(I);
v=double(I);
vec=[1:256];
l=0;
for k=0:255 
    for i=1:nl
        for j=1:nc
            if v(i,j)==k 
               l=l+1;
            end
        end
    end
    vec(k+1)=l;
    l=0;
end
axes(handles.imgT);subimage(plot(vec));



% --------------------------------------------------------------------
function Moyenneur33_Callback(hObject, eventdata, handles)
% hObject    handle to Moyenneur33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[n,m]=size(image);
image = double(image);
b=image;
H=(1/9)*[1 1 1 ; 1 1 1 ; 1 1 1 ];
for x = 2 : n-1
    for y = 2 : m-1
     %img(x,y)=mean([image(x-1,y-1)+image(x-1,y)+image(x-1,y+1)+image(x,y-1)+image(x,y)+image(x,y+1)+image(x+1,y-1)+image(x+1,y)+image(x+1,y+1)]);
      f=image(x-1:x+1,y-1:y+1);
      v=f.*H;
      %b=conv2(img,H);
      b(x,y)=sum(v(:));
      %b(x,y)=mean(f(:));
    end 
end
b=uint8(b);
%imshow(b);
axes(handles.imgT);
subimage(b);
title('Image filtrée Moy3');


% --------------------------------------------------------------------
function Moyenneur55_Callback(hObject, eventdata, handles)
% hObject    handle to Moyenneur55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[n,m]=size(image);
image = double(image);
b=image;
H=(1/25)*[1 1 1 1 1 ; 1 1 1 1 1 ; 1 1 1 1 1 ; 1 1 1 1 1 ; 1 1 1 1 1];
for x = 3 : n-2
    for y = 3 : m-2
     f=image(x-2:x+2,y-2:y+2);
      v=f.*H;
      b(x,y)=sum(v(:));
    end 
end
b=uint8(b);
axes(handles.imgT);
     subimage(b);
handles.ima_traite = b;
handles.output = hObject;
guidata(hObject, handles);
title('Image filtrée Moy5');


% --------------------------------------------------------------------
function Gaussien55_Callback(hObject, eventdata, handles)
% hObject    handle to Gaussien55 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[n,m]=size(image);
image = double(image);
b=image;
H=(1/256)*[1 4 6 4 1 ; 4 16 24 16 4 ; 6 24 36 24 6 ; 4 16 24 16 4 ; 1 4 6 4 1];
for x = 3 : n-2
    for y = 3 : m-2
  f=image(x-2:x+2,y-2:y+2);
      v=f.*H;
      b(x,y)=sum(v(:));
    end 
end
b=uint8(b);
axes(handles.imgT);
subimage(b);
handles.ima_traite = b;
handles.output = hObject;
guidata(hObject, handles);
title('Image filtrée Gauss5');



% --------------------------------------------------------------------
function Gaussien33_Callback(hObject, eventdata, handles)
% hObject    handle to Gaussien33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[n,m]=size(image);
image = double(image);
b=image;
H=(1/16)*[1 2 1 ;2 4 2 ; 1 2 1];
for x = 2 : n-1
    for y = 2 : m-1
    f=image(x-1:x+1,y-1:y+1);
      v=f.*H;
      b(x,y)=sum(v(:));
    end 
end
b=uint8(b);
axes(handles.imgT);
subimage(b);
handles.ima_traite = b;
handles.output = hObject;
guidata(hObject, handles);
title('Image filtrée Gauss3');


% --------------------------------------------------------------------
function Contraste_Callback(hObject, eventdata, handles)
% hObject    handle to Contraste (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
%[n,m]=size(image);
image = double(image);
%output=image;

%ima=imread('cameraman.tif');
[l c]=size(image);
image = double(image);
v=image;
for i=1:l
    for j=1:c
      fpixel = (image(i,j)-128)*5 + 128; 
    % on vérifie que la valeur obtenue est bien dans [0..255]
    if( fpixel>255 )
      fpixel = 255;
    else if( fpixel<0 )
      fpixel = 0;
        end 
    end
    
   v(i,j) = fpixel;
    end
end  
v=uint8(v); 
axes(handles.imgT);
subimage(v);
title('Image contraste');


% --------------------------------------------------------------------
function Binarisation_Callback(hObject, eventdata, handles)
% hObject    handle to Binarisation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I4 = handles.courant_data;
%****************************
% Calcule du seuil
%****************************
%calucle des m:
m0=1;
m1=mean2(I4);
m2=mean2(I4.^2);
m3=mean2(I4.^3);
%calcule des C:
C1=(m3-(m1*m2))/(m2-m1);
C0=(-m2-(C1*m1))/m0;
%calcule des z:
z1=(-C1-sqrt(C1^2-4*C0))/2;
z2=(-C1+sqrt(C1^2-4*C0))/2;
seuil=(z1+z2)/2;
% première solution:
% bin=zeros(240,320);
% for i=1:240
% for j=1:320
% if I4(i,j)>seuil;
% bin(i,j)=255;
% end
% end
% end
bin=(I4>seuil)*255;
text(2,10,num2str(seuil));
%
% Solution via matlab de la binarisation :
% level=graythresh(I4); %calcule seuil
% bin = im2bw(I4,level); %binarisation matlab
axes(handles.imgO);
subimage(I4);
axes(handles.imgT);
handles.ima_traite = bin;
subimage(handles.ima_traite);

%Grrr
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
title('Image binaire');


% --------------------------------------------------------------------
function Laplacien_Callback(hObject, eventdata, handles)
% hObject    handle to Laplacien (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;

% Convertir l'image en niveaux de gris
image = rgb2gray(image);

%image=imnoise(imageO,'salt & pepper', 0.05);
[n,m]=size(image);
image = double(image);
%b=image;
[n m]=size(image);
b=zeros(n,m); %initialize l 'img de sortie par Zero car la resulta est binaire
%M1=[0 1 0;1 -4 1;0 1 0];
M1=[-1 -1 -1;-1 8 -1;-1 -1 -1];
for i=2:n-1
    for j=2:m-1
        V=image((i-1:i+1),(j-1:j+1));
        S=V.*M1;
        b(i,j)=sum(S(:));
    end
end
b=uint8(b);
axes(handles.imgT);
     subimage(b);
     title('Image filtrée Laplacien');


% --------------------------------------------------------------------
function Gradient_Callback(hObject, eventdata, handles)
% hObject    handle to Gradient (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[n,m]=size(image);
image = double(image);
output=image;% inistalize l'entrer par la sortie pour éviter le probléme  des bornes
%image=rgb2gray(image);

[m,n] = size(image);
output=zeros(size(image)); 
outputhor=zeros(size(image)); 
outputver=zeros(size(image)); 
maskhor = [0,0,0;-1,0,1;0,0,0]; %maske horizentale
maskver = [0,-1,0;0,0,0;0,1,0]; % maske verticale 
for i=4:(m-3)
   for j=4:(n-3) 
      for k=1:3         
          for l=1:3
            outputhor(i,j) = outputhor(i,j)+image(i-k,j-l)*maskhor(k,l);            
            outputver(i,j) = outputver(i,j)+image(i-k,j-l)*maskver(k,l);          
          end
      end
    end
 end
%mymin=min(min(output))
%mymax=max(max(output))
for i=4:(m-3)
for j=4:(n-3)       
    output(i,j)=sqrt(outputhor(i,j)*outputhor(i,j) + outputver(i,j)*outputver(i,j));
end 
end 
%outputhor=uint8(outputhor); 
%outputver=uint8(outputver); 
output=uint8(output); 

%b=uint8(b);
axes(handles.imgT);
subimage(output);
title('Image filtrée Grad');

%figure(10);colormap(gray(256));imagesc(outputhor);title('gradient hor'); 
%figure(11);colormap(gray(256));imagesc(outputver);title('gradient ver'); 
%figure(12);colormap(gray(256));imagesc(output);title('gradient');


% --------------------------------------------------------------------
function Sobel_Callback(hObject, eventdata, handles)
% hObject    handle to Sobel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;

% Convertir l'image en niveaux de gris
image = rgb2gray(image);

[n,m]=size(image);
image = double(image);
output=image;
%image=rgb2gray(image);

output=zeros(size(image)); 
outputhor=zeros(size(image)); 
outputver=zeros(size(image)); 
%maskhor = [-1,0,1;-b,0,b;-1,0,1]; 
%maskver = [-1,-b,-1;0,0,0;1,b,1];
%b=2 pour sobel
maskhor = [-1,0,1;-2,0,2;-1,0,1]; 
maskver = [-1,-2,-1;0,0,0;1,2,1];

for i=4:(m-3)
   for j=4:(n-3) 
      for k=1:3          
          for l=1:3
            outputhor(i,j) = outputhor(i,j)+image(i-k,j-l)*maskhor(k,l);             
            outputver(i,j) = outputver(i,j)+image(i-k,j-l)*maskver(k,l);          
          end
      end
    end
 end
%mymin=min(min(output))
%mymax=max(max(output))
for i=4:(m-3)
   for j=4:(n-3) 
output(i,j)=sqrt(outputhor(i,j)*outputhor(i,j) + outputver(i,j)*outputver(i,j)); 
   end
end
%outputhor=uint8(outputhor); 
%outputver=uint8(outputver); 
output=uint8(output); 
axes(handles.imgT);
subimage(output);
title('Image filtrée sobel');

%figure(10);colormap(gray(256));imagesc(outputhor);title(gradient hor); 
%figure(11);colormap(gray(256));imagesc(output);title(gradient ver); 
%figure(12);colormap(gray(256));imagesc(output);title(gradient);


% --------------------------------------------------------------------
function Prewitt_Callback(hObject, eventdata, handles)
% hObject    handle to Prewitt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;

% Convertir l'image en niveaux de gris
image = rgb2gray(image);

[n,m]=size(image);
image = double(image);
output=image;
%image=rgb2gray(image);

output=zeros(size(image)); 
outputhor=zeros(size(image)); 
outputver=zeros(size(image)); 
%maskhor = [-1,0,1;-b,0,b;-1,0,1]; 
%maskver = [-1,-b,-1;0,0,0;1,b,1];
%b=1 pour prewitt
maskhor = [-1,0,1;-1,0,1;-1,0,1]; 
maskver = [-1,-1,-1;0,0,0;1,1,1];

for i=4:(m-3)
   for j=4:(n-3) 
      for k=1:3          
          for l=1:3
            outputhor(i,j) = outputhor(i,j)+image(i-k,j-l)*maskhor(k,l);             
            outputver(i,j) = outputver(i,j)+image(i-k,j-l)*maskver(k,l);          
          end
      end
    end
 end
%mymin=min(min(output))
%mymax=max(max(output))
for i=4:(m-3)
   for j=4:(n-3) 
output(i,j)=sqrt(outputhor(i,j)*outputhor(i,j) + outputver(i,j)*outputver(i,j)); 
   end
end
%outputhor=uint8(outputhor); 
%outputver=uint8(outputver); 
output=uint8(output); 
axes(handles.imgT);
subimage(output);
title('Image filtrée Prew');

%figure(10);colormap(gray(256));imagesc(outputhor);title(gradient hor); 
%figure(11);colormap(gray(256));imagesc(output);title(gradient ver); 
%figure(12);colormap(gray(256));imagesc(output);title(gradient);



% --------------------------------------------------------------------
function Robert_Callback(hObject, eventdata, handles)
% hObject    handle to Robert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = handles.courant_data;

% Convertir l'image en niveaux de gris
image = rgb2gray(image);

[n,m]=size(image);
image = double(image);
 %num = get(handles.slider1, 'value');
% set(handles.edit1, 'String', num2str(num));
for x=1:n-1
 for y=1:m-1
  b(x,y)= abs(uint8( double(image(x,y))-double(image(x+1,y+1))))+ abs(uint8( double(image(x,y+1)) - double(image(x+1,y))));
 end
end
    % num = get(handles.slider1, 'Value');
    % set(handles.txt1, 'String', num2str(num));
        %Seuillage
        [n,m]=size(image);
        for i=1:n-1
         for j=1:m-1
          if b(i,j) < 25
            b(i,j)=0;
          end
         end
        end
           %
  handles.ima_traite = b;
  axes(handles.imgT);
  subimage(b);
%Grrr
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
title('Image filtrée Rob');


% --------------------------------------------------------------------
function FPB_Callback(hObject, eventdata, handles)
% hObject    handle to FPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

I = handles.courant_data;

 F=fftshift(fft2(I)); 
% %calcul de la taille de l'image; 
 M=size(F,1); 
 N=size(F,2); 
 P=size(F,3); 
 H0=zeros(M,N); 
 D0=1; 
 M2=round(M/2); 
 N2=round(N/2); 
 H0(M2-D0:M2+D0,N2-D0:N2+D0)=1; 
 for i=1:M 
  for j=1:N 
      G(i,j)=F(i,j)*H0(i,j); 
 end 
 end 
 g=ifft2(G); 
 imshow(abs(g),[0,255]);
 title('Image filtrée Fpb');


% --------------------------------------------------------------------
function FPH_Callback(hObject, eventdata, handles)
% hObject    handle to FPH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I=handles.courant_data;
%charge; 
F=fftshift(fft2(I)); 
%calcul de la taille de l'image; 
M=size(F,1); 
N=size(F,2); 
P=size(F,3); 

H1=ones(M,N); 
D0=1; 
M2=round(M/2); 
N2=round(N/2); 
H1(M2-D0:M2+D0,N2-D0:N2+D0)=0; 
for i=1:M 
for j=1:N 
G(i,j)=F(i,j)*H1(i,j); 
end 
end 
g=ifft2(G); 
%subplot(1,2,1);imshow(I);title('image originale'); 
%subplot(1,2,2);
imshow(255-abs(g),[0,255]);
title('Image filtrée Fph');
%title('image filtrée');


% --------------------------------------------------------------------
function Conique_Callback(hObject, eventdata, handles)
% hObject    handle to Conique (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[ n,m]=size(image);
image = double(image);
b=image;
H=(1/25)*[0 0 1 0 0 ; 0 2 2 2 0 ; 1 2 5 2 1 ; 0 2 2 2 0 ; 0 0 1 0 0];
for x = 3 : n-2
    for y = 3 : m-2
       f=image(x-2:x+2,y-2:y+2);
      v=f.*H;
      b(x,y)=sum(v(:));
    end 
end
b=uint8(b);
axes(handles.imgT);
subimage(b);
handles.ima_traite = b;
handles.output = hObject;
guidata(hObject, handles);
title('Image filtrée Conique');


% --------------------------------------------------------------------
function Pyramidal_Callback(hObject, eventdata, handles)
% hObject    handle to Pyramidal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image=handles.courant_data;
[n,m]=size(image);
image = double(image);
b=image;
H=(1/81)*[1 2 3 2 1 ; 2 4 6 4 2 ; 3 6 9 6 3 ; 2 4 6 4 2 ; 1 2 3 2 1];
for x = 3 : n-2
    for y = 3 : m-2
          f=image(x-2:x+2,y-2:y+2);
      v=f.*H;
      b(x,y)=sum(v(:));
    end 
end
b=uint8(b);
axes(handles.imgT);
subimage(b);
handles.ima_traite = b;
handles.output = hObject;
guidata(hObject, handles);
title('Image filtrée Payr');


% --------------------------------------------------------------------
function Histogramme_Callback(hObject, eventdata, handles)
% hObject    handle to Histogramme (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA) 
 img =handles.courant_data;
      d = length(size(img));
      if d==3
          I=rgb2gray(img);
      elseif d==2 
          I = img;
      end
 axes(handles.imgO);
 imshow(I);
 
 [nl,nc]=size(I);
 v=double(I);
 vec=[1:256];
 l=0;
 for k=0:255
     for i=1:nl
         for j=1:nc
             if v(i,j)==k
                 l=l+1;
             end
         end
     end
     vec(k+1)=l;
     l=0;
 end
 axes(handles.imgT);plot(vec);
 title('histogramme de L Image');
 


% --------------------------------------------------------------------
function Morphologie_math_Callback(hObject, eventdata, handles)
% hObject    handle to Morphologie_math (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Erosion_Callback(hObject, eventdata, handles)
% hObject    handle to Erosion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% read image
I=handles.courant_data;
%I=imread('lenna.png');  
 
% convert to binary 
I=im2bw(I);
 
% create structuring element             
se=ones(5, 5);
 
% store number of rows
% in P and number of columns in Q.           
[P, Q]=size(se);
 
% create a zero matrix of size I.       
In=zeros(size(I, 1), size(I, 2));
 
for i=ceil(P/2):size(I, 1)-floor(P/2)
    for j=ceil(Q/2):size(I, 2)-floor(Q/2)
 
        % take all the neighbourhoods.
        on=I(i-floor(P/2):i+floor(P/2), j-floor(Q/2):j+floor(Q/2));
 
        % take logical se
        nh=on(logical(se));
       
        % compare and take minimum value of the neighbor
        % and set the pixel value to that minimum value.
        In(i, j)=min(nh(:));     
    end
end
 
imshow(In);
title('Image érode');

% --------------------------------------------------------------------
function Dilatation_Callback(hObject, eventdata, handles)
% hObject    handle to Dilatation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% read image
I=handles.courant_data;
%I=imread('lenna.png');    
 
% convert to binary
I=im2bw(I); 
 
% create structuring element            
se=ones(5, 5);
 
% store number of rows in P and
% number of columns in Q.          
[P, Q]=size(se);
 
% create a zero matrix of size I.       
In=zeros(size(I, 1), size(I, 2));
 
for i=ceil(P/2):size(I, 1)-floor(P/2)
    for j=ceil(Q/2):size(I, 2)-floor(Q/2)
 
        % take all the neighbourhoods.
        on=I(i-floor(P/2):i+floor(P/2), j-floor(Q/2):j+floor(Q/2)); 
        
        % take logical se
        nh=on(logical(se));   
 
        % compare and take minimum value of the neighbor
        % and set the pixel value to that minimum value.   
        In(i, j)=max(nh(:));     
    end
end
 
imshow(In);
title('Image dilaté');

% --------------------------------------------------------------------
function Ouverture_open_Callback(hObject, eventdata, handles)
% hObject    handle to Ouverture_open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function Fermeture_Close_Callback(hObject, eventdata, handles)
% hObject    handle to Fermeture_Close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function Open_pre_Callback(hObject, eventdata, handles)
% hObject    handle to Open_pre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

I=handles.courant_data;   
 
% convert to binary
I=im2bw(I);
% create structuring element
se = strel('square', 3);

% Perform opening on the binary image
opened_image = imopen(I, se);
imshow(opened_image);
title('Image opened');


% --------------------------------------------------------------------
function Close_pre_Callback(hObject, eventdata, handles)
% hObject    handle to Close_pre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I=handles.courant_data;
%I=imread('lenna.png');    
 
% convert to binary
I=im2bw(I);
% Define the structuring element for closing
se = strel('square', 3);

% Perform opening on the binary image
closed_image = imclose(I, se);
imshow(closed_image);
title('Image closed');


% --------------------------------------------------------------------
function Gradient_interne_Callback(hObject, eventdata, handles)
% hObject    handle to Gradient_interne (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

I=handles.courant_data;
%I=imread('lenna.png');    
 
% convert to binary
I=im2bw(I);
% Define the structuring element for closing
se = strel('square', 3);

% Perform opening on the binary image
contour_interne = I - imerode(I, se);
imshow(contour_interne);
title('contour interne');

% --------------------------------------------------------------------
function Gradient_externe_Callback(hObject, eventdata, handles)
% hObject    handle to Gradient_externe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I=handles.courant_data;
%I=imread('lenna.png');    
 
% convert to binary
I=im2bw(I);
% Define the structuring element for closing
se = strel('square', 3);

% Perform opening on the binary image
contour_externe =imdilate(I, se)-I;
imshow(contour_externe);
title('contour externe');

% --------------------------------------------------------------------
function Gradient_contour_Callback(hObject, eventdata, handles)
% hObject    handle to Gradient_contour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Gradient_morphologique_Callback(hObject, eventdata, handles)
% hObject    handle to Gradient_morphologique (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I=handles.courant_data;
%I=imread('lenna.png');    
 
% convert to binary
I=im2bw(I);
% Define the structuring element for closing
se = strel('square', 3);

% Perform opening on the binary image
contour_morph =imdilate(I, se)-imerode(I, se);
imshow(contour_morph);
title('contour morphologique');

% --------------------------------------------------------------------
function Untitled_43_Callback(hObject, eventdata, handles)
% hObject    handle to Close_pre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Ouverture_Callback(hObject, eventdata, handles)
% hObject    handle to Ouverture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I=handles.courant_data;    
%se = strel('square', 3);
% create structuring element    
%se = strel('square', 3);
se=ones(5, 5);
% convert to binary
I=im2bw(I);

% Perform morphological opening
img_opened = imerode(imdilate(I, se), se);
imshow(img_opened);
title('image opened');


% --------------------------------------------------------------------
function Fermeture_Callback(hObject, eventdata, handles)
% hObject    handle to Fermeture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I=handles.courant_data;
% create structuring element    
%se = strel('square', 3);
se=ones(5, 5);
% convert to binary
I=im2bw(I);
% Perform morphological opening
img_closed = imdilate(imerode(I, se), se);
imshow(img_closed);
title('image closed');


% --------------------------------------------------------------------
function Points_interet_Callback(hObject, eventdata, handles)
% hObject    handle to Points_interet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Det_SUSAN_Callback(hObject, eventdata, handles)
% hObject    handle to Det_SUSAN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Load the image and convert it to grayscale
img=handles.courant_data;
img = rgb2gray(img);

% Apply Gaussian smoothing filter
h = fspecial('gaussian', [5 5], 2);
img_smoothed = imfilter(img, h);

% Compute gradient magnitude
[Gx,Gy] = imgradientxy(img_smoothed, 'sobel');
grad_mag = sqrt(Gx.^2 + Gy.^2);

% Initialize the edge map
edge_map = zeros(size(img));

% Loop through each pixel
for i = 1:size(img,1)
    for j = 1:size(img,2)
        % Get the intensity value of the current pixel
        p = img_smoothed(i,j);
        
        % Get the pixels in a circular region around the current pixel
        circle = getCircularRegion(i,j,img_smoothed,20);
        
        % Count the number of pixels with intensity values close to p
        count = sum(abs(circle - p) < 30);
        
        % If the count is below a threshold, mark the pixel as an edge
        if count < 20
            edge_map(i,j) = 1;
        end
    end
end

% Display the edge map

plot(edge_map);
% --------------------------------------------------------------------
function HARRIS_Callback(hObject, eventdata, handles)
% hObject    handle to HARRIS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img=handles.courant_data;

%img = rgb2gray(img);
imd = double(img);
[m, n]=size(imd)%cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
sigma=1; k=0.04; w=5*sigma; seuil=50; r=6; sze=2*r+1;
dx=[-1 0 1]; dy=dx'; % filtre dérivatif
g = fspecial('gaussian',max(1,fix(w)), sigma); % filtre gaussien
Ix=conv2(imd,dx,'same');
Iy=conv2(imd,dy,'same');
Ix2 = conv2(Ix.^2, g, 'same'); 
Iy2 = conv2(Iy.^2, g, 'same');
Ixy = conv2(Ix.*Iy, g, 'same');
R=Ix2.*Iy2-Ixy.^2-k*(Ix2+Iy2).^2;
R1=(1000/(max(max(R))))*R;
%****** Seuillage et extraction des points d'intérêt ********
[u,v]=find(R1<=seuil);
nb=length(u);
for l=1:nb
R1(u(l),v(l))=0;
end
R11=zeros(m+2*r,n+2*r);
R11(r+1:m+r,r+1:n+r)=R1;
[m1,n1]=size(R11);
for i=r+1:m1-r
for j=r+1:n1-r
fenetre=R11(i-r:i+r,j-r:j+r);
ma=max(max(fenetre));
if fenetre(r+1,r+1)<ma
R11(i,j)=0;
end
end
end
R11=R11(r+1:m+r,r+1:n+r);
[x,y]=find(R11);
%===============================================================
MX=ordfilt2(R1,sze^2,ones(sze));
R11 = (R1==MX)&(R1>seuil);
[x,y]=find(R11);
%************* affichage des points d'intérêt **************
nb=length(x);
imshow(img)
hold on
plot(y,x,'r.')

% --------------------------------------------------------------------
function Points_d_interet_Callback(hObject, eventdata, handles)
% hObject    handle to Points_d_interet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function SUSAN_Callback(hObject, eventdata, handles)
% hObject    handle to SUSAN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
im=handles.courant_data;
%im = rgb2gray(im);
d = length(size(im));

% =======================conversion de l'image=============================

    if d==3
        image=double(rgb2gray(im));
    elseif d==2
        image=double(im);
    end
[n,m]=size(image);
% =============================données=====================================
rayon=1;
alpha=80;
r=5;
alpha=alpha/100;
% ========================génerateur de mask=============================
mask=zeros(2*rayon+1);
b=ones(rayon+1);
for i=1:rayon+1
    for j=1:rayon+1
        if (rayon==1)
            if(j>i)
            b(i,j)=0;
            end
        else
            if(j>i+1)
            b(i,j)=0;
            end
        end
    end
end
mask(1:rayon+1,rayon+1:2*rayon+1)=b;
mask(1:rayon+1,1:rayon+1)=rot90(b);
mask0=mask;
mask0=flip(mask0,1);
mask=mask0+mask;
mask(rayon+1,:)=mask(rayon+1,:)-1;
% ==========================réponse maximale============================
max_reponse=sum(sum(mask));
% =====================balayage de toute l'image===========================
f=zeros(n,m);
for i=(rayon+1):n-rayon
    for j=(rayon+1):m-rayon
        image_courant=image(i-rayon:i+rayon,j-rayon:j+rayon);
        image_courant_mask=image_courant.*mask;
        inteniste_cental= image_courant_mask(rayon+1,rayon+1);
        s=exp(-1*(((image_courant_mask-inteniste_cental)/max_reponse).^6));
        somme=sum(sum(s));
% si le centre du mask est un 0 il faut soustraire les zeros des filtres
       if (inteniste_cental==0)
        somme=somme-length((find(mask==0)));
       end
        f(i,j)=somme;
    end
end

% =============selection et seuillage des points d'interét=================
ff=f(rayon+1:n-(rayon+1),rayon+1:m-(rayon+1));
minf=min(min(ff));
maxf=max(max(f));
fff=f;
d=2*r+1;
temp1=round(n/d);
   if (temp1-n/d)<0.5 &&(temp1-n/d)>0
    temp1=temp1-1;
   end
    temp2=round(m/d);
   if (temp2-m/d)<0.5 &&(temp2-m/d)>0
    temp2=temp2-1;
   end
    fff(n:temp1*d+d,m:temp2*d+d)=0;

for i=(r+1):d:temp1*d+d
   for j=(r+1):d:temp2*d+d
        window=fff(i-r:i+r,j-r:j+r);
        window0=window;
        [xx,yy]=find(window0==0);
        for k=1:length(xx)
            window0(xx(k),yy(k))=max(max(window0));
        end
        minwindow=min(min(window0));
        [y,x]=find(minwindow~=window & window<=minf+alpha*(maxf-minf) & window>0);
        [u,v]=find(minwindow==window);
    if length(u)>1
        for l=2:length(u)
            fff(i-r-1+u(l),j-r-1+v(l))=0 ;
        end
    end
   if ~isempty(x)
    for l=1:length(y)
        fff(i-r-1+y(l),j-r-1+x(l))=0 ;
    end
   end
   end
end
seuil=minf+alpha*(maxf-minf);
[u,v]=find(minf<=fff & fff<=seuil );

% ==============affichage des resultats====================================

%axes(handles.imgS);
imshow(im);
hold on;
plot(v,u,'.r','MarkerSize',10);
hold off;
message = sprintf(' le nombre des points d''intérêts: %d      ',length(v));
msgbox(message);


% --------------------------------------------------------------------
function S_Callback(hObject, eventdata, handles)
% hObject    handle to S (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function SS_Callback(hObject, eventdata, handles)
% hObject    handle to SS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Load an image and convert to binary

I = handles.courant_data;
I = rgb2gray(I);
BW = edge(I, 'canny');

% Initialize accumulator array
d = sqrt(size(I,1)^2 + size(I,2)^2);
nbins_r = round(d);
nbins_theta = 180;
accumulator = zeros(nbins_r, nbins_theta);

% Define circle voting function
[X, Y] = meshgrid(1:size(I,2), 1:size(I,1));
vote_circle = @(x, y, r, theta) round(x + r * cosd(theta)) + size(I,2) * (round(y + r * sind(theta)) - 1);

% Hough transform
for x = 1:size(BW,2)
    for y = 1:size(BW,1)
        if BW(y,x)
            for theta = 1:nbins_theta
                r = round(x * cosd(theta) + y * sind(theta));
                accumulator(r, theta) = accumulator(r, theta) + 1;
            end
        end
    end
end

% Identify local maxima
local_max = imregionalmax(accumulator);
[r_idx, theta_idx] = find(local_max);

% Overlay circles on original image
imshow(I); hold on;
for i = 1:length(r_idx)
    x = r_idx(i) * cosd(theta_idx(i)) - size(I,2) / 2;
    y = r_idx(i) * sind(theta_idx(i)) - size(I,1) / 2;
    r = sqrt(x^2 + y^2);
    draw_circle(x, y, r);
end
hold off;

function draw_circle(x, y, r)
    t = linspace(0, 2*pi, 100);
    X = x + r * cos(t);
    Y = y + r * sin(t);
    plot(X, Y, 'r', 'LineWidth', 2);



% --------------------------------------------------------------------
function HH_Callback(hObject, eventdata, handles)
% hObject    handle to HH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function SUS_Callback(hObject, eventdata, handles)
% hObject    handle to SUS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Load the image
I = handles.courant_data;
% Load the image
%I = imread('image.jpg');
I = rgb2gray(I);

% Define the SUSAN threshold
threshold = 50;

% Initialize the output image
output = zeros(size(I));

% Loop over each pixel in the image
for i = 3:size(I,1)-2
    for j = 3:size(I,2)-2
        % Extract the 37x37 neighborhood around the current pixel
        neighborhood = I(i-18:i+18, j-18:j+18);
        
        % Calculate the SUSAN function for the current pixel
        g = sum(abs(neighborhood(:) - I(i,j)) < threshold);
        
        % Threshold the SUSAN function
        if g < 20
            output(i,j) = 0;
        else
            output(i,j) = 255;
        end
    end
end

% Display the output image
imshow(output, [0 255]);


% --------------------------------------------------------------------
function Modele_electrostatique_Callback(hObject, eventdata, handles)
% hObject    handle to Modele_electrostatique (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img=handles.courant_data;
k=0.04; sigma=1; seuil=100; r=6; w=5*sigma;
[m,n]=size(img); imd=double(img);
dxa=[-sqrt(2)/4 0 sqrt(2)/4 ; -1 0 1 ; -sqrt(2)/4 0 sqrt(2)/4];
% dxa=[sqrt(2)/4 0 -sqrt(2)/4; 1 0 -1; sqrt(2)/4 0 -sqrt(2)/4];
dya=dxa'; % derivée verticale
g=fspecial('gaussian',max(1,fix(5*sigma)),sigma); % gaussien
Ixa=conv2(imd,dxa,'same');
Iya=conv2(imd,dya,'same');
Ixa2 = conv2(Ixa.^2, g, 'same'); 
Iya2 = conv2(Iya.^2, g, 'same');
Ixya = conv2(Ixa.*Iya, g,'same');
R=Ixa2.*Iya2-Ixya.^2-k*(Ixa2+Iya2).^2;
R1=(1000/(max(max(R))))*R; %normalisation
[u,v]=find(R1<=seuil);
nb=length(u);
for k=1:nb
R1(u(k),v(k))=0;
end
R11=zeros(m+2*r,n+2*r);
R11(r+1:m+r,r+1:n+r)=R1;
[m1,n1]=size(R11);
for i=r+1:m1-r
for j=r+1:n1-r
fenetre=R11(i-r:i+r,j-r:j+r);
ma=max(max(fenetre));
if fenetre(r+1,r+1)<ma
R11(i,j)=0;
end
end
end
%subplot(1,2,2);cccccccccccccccccccccccccccccccc
imshow(img)
hold on
R11=R11(r+1:m+r,r+1:n+r);
[x,y]=find(R11);
nb=length(x)
plot(y,x,'.r');


% --------------------------------------------------------------------
function fft_Callback(hObject, eventdata, handles)
% hObject    handle to fft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function fft_auto_Callback(hObject, eventdata, handles)
% hObject    handle to fft_auto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Load the image
img = handles.courant_data;
%img = rgb2gray(img);

% Calculate the Fourier transform
F = fft2(img);

% Compute magnitude and phase of the Fourier transform
F_magnitude = abs(F);

FF=log(1 + F_magnitude)
% Plot the original image and the Fourier transform
imshow(FF, []);


% --------------------------------------------------------------------
function Filtrage_Morphologique_Callback(hObject, eventdata, handles)
% hObject    handle to Filtrage_Morphologique (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Lire l'image
img = handles.courant_data;

% Convertir l'image en niveaux de gris
img_gray = rgb2gray(img);
%img_gray=im2bw(img);
% Définir la taille de l'élément structurant
strel_size = 3;
%se = strel('disk', strel_size);
se = strel('square', strel_size);

% Appliquer une opération d'ouverture à l'image
img_opened = imopen(img_gray, se);

% Appliquer une opération de fermeture à l'image
img_closed = imclose(img_opened, se);
% Afficher l'image lissée
imshow(img_closed);
% --------------------------------------------------------------------
function Detection_Pics_Callback(hObject, eventdata, handles)
% hObject    handle to Detection_Pics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Egalisation_Callback(hObject, eventdata, handles)
% hObject    handle to Egalisation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = handles.courant_data;
numPixels = numel(image); %Calculez le nombre de pixels dans l'image en utilisant la fonction "numel" :
histogram = zeros(1,256);
%Parcourez chaque pixel de l'image et incrémentez l'élément correspondant dans le vecteur d'histogramme :
for i = 1:numPixels
    pixelValue = image(i);
    histogram(pixelValue+1) = histogram(pixelValue+1) + 1;
end
cdf = cumsum(histogram) / numPixels; %Calculez la fonction de distribution cumulée (CDF) du vecteur d'histogramme :
cdf = cdf * 255; % Multipliez la CDF par le nombre maximum de niveaux de gris (255) pour normaliser la CDF :
cdf = uint8(cdf);%Convertissez la CDF en entiers non signés de 8 bits :
outputImage = cdf(image+1); %Appliquez la transformation de l'histogramme en utilisant la CDF :
imshow(outputImage);


% --------------------------------------------------------------------
function Kirsch_Callback(hObject, eventdata, handles)
% hObject    handle to Kirsch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Charger l'image
img =  handles.courant_data;

% Convertir l'image en niveaux de gris
img_gray = rgb2gray(img);

% Définir le noyau de Kirsh (matrice de convolution)
kernel = [-3 -3 5; -3 0 5; -3 -3 5];

% Appliquer la convolution à l'image avec le noyau de Kirsh
filtered_img = conv2(double(img_gray), double(kernel), 'same');

% Afficher l'image filtrée
imshow(filtered_img, []), title('image par le filtre de Kirsch');
% --------------------------------------------------------------------
function Marr_hildreth_Callback(hObject, eventdata, handles)
% hObject    handle to Marr_hildreth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Lecture de l'image et conversion en niveaux de gris
img = handles.courant_data;
img = rgb2gray(img);

% Définition des paramètres du filtre
sigma = 2; % écart type du filtre gaussien
threshold = 0.1; % seuil pour la détection des contours

% Calcul du filtre de Laplacien de Gaussien
hsize = ceil(sigma*3)*2+1; % taille du noyau gaussien
gauss = fspecial('gaussian', hsize, sigma);
lap = fspecial('laplacian', 0);
filtre = conv2(gauss, lap, 'same');

% Filtrage de l'image avec le filtre de Laplacien de Gaussien
img_filtree = conv2(double(img), filtre, 'same');

% Recherche des passages à zéro du filtre
zero_crossing = zeros(size(img_filtree));
for i=2:size(img_filtree,1)-1
    for j=2:size(img_filtree,2)-1
        neighbours = img_filtree(i-1:i+1,j-1:j+1);
        if img_filtree(i,j) > threshold && any(neighbours(:) < -threshold)
            zero_crossing(i,j) = 1;
        end
    end
end

% Affichage des contours détectés

imshow(zero_crossing), title('Contours détectés');

% --------------------------------------------------------------------
function Sombre_Callback(hObject, eventdata, handles)
% hObject    handle to Sombre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Charger l'image
img = handles.courant_data;

% Convertir l'image en niveaux de gris
img_gray = rgb2gray(img);

% Inverser l'image pour obtenir les pics sombres
inverted_img = imcomplement(img_gray);

% Créer un élément structurant pour la détection des pics sombres
se = strel('disk', 5);

% Appliquer une opération d'érosion à l'image inversée
eroded_img = imerode(inverted_img, se);

% Appliquer une opération de dilatation à l'image érodée
dilated_img = imdilate(eroded_img, se);

% Soustraire l'image érodée de l'image originale pour obtenir l'image des pics sombres
dark_peak_img = img_gray - eroded_img;

% Afficher l'image avec les pics sombres détectés
imshow(dark_peak_img);
title('image Detection pics Sombre');
% --------------------------------------------------------------------
function claire_Callback(hObject, eventdata, handles)
% hObject    handle to claire (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Charger l'image
img =  handles.courant_data;

% Convertir l'image en niveaux de gris
img_gray = rgb2gray(img);

% Créer un élément structurant pour la détection des pics clairs
se = strel('disk', 5);

% Appliquer une opération d'érosion à l'image originale
eroded_img = imerode(img_gray, se);

% Appliquer une opération de dilatation à l'image érodée
dilated_img = imdilate(eroded_img, se);

% Soustraire l'image érodée de l'image originale pour obtenir l'image des pics clairs
bright_peak_img = img_gray - eroded_img;

% Afficher l'image avec les pics clairs détectés
imshow(bright_peak_img);
title('image Detection Pics Claire');


% --------------------------------------------------------------------


% --------------------------------------------------------------------
function Canny_Callback(hObject, eventdata, handles)
% hObject    handle to Canny (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
img = handles.courant_data;
% Convert the image to grayscale
img_gray = rgb2gray(img);

% Apply the Canny edge detection algorithm with default parameters
%edges = edge(img_gray, 'canny');
low_threshold = 0.1;
high_threshold = 0.3;
edges = edge(img_gray, 'canny', [low_threshold high_threshold]);

% Display the result
imshow(edges);
title('image filtre Canny');


% --------------------------------------------------------------------
function Top_Hat_Callback(hObject, eventdata, handles)
% hObject    handle to Top_Hat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function White_Callback(hObject, eventdata, handles)
% hObject    handle to White (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Charger l'image
img = handles.courant_data;

% Convertir l'image en niveaux de gris
img_gray = rgb2gray(img);

% Définir un élément structurant pour l'ouverture morphologique
se = strel('disk', 5);

% Effectuer l'ouverture morphologique
img_opened = imopen(img_gray, se);

% Soustraire l'image originale de l'image ouverte
img_tophat = img_gray - img_opened;
% Display the result
imshow(img_tophat);
title('image WHITE TOP HAT');



% --------------------------------------------------------------------
function Black_Callback(hObject, eventdata, handles)
% hObject    handle to Black (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Charger l'image
img = handles.courant_data;

% Convertir l'image en niveaux de gris
img_gray = rgb2gray(img);

% Définir un élément structurant pour l'ouverture morphologique
se = strel('disk', 5);

% Effectuer l'ouverture morphologique
img_closed = imclose(img_gray, se);

% Soustraire l'image originale de l'image ouverte
img_tophat = img_closed - img_gray;
% Display the result
imshow(img_tophat);
title('image Black TOP HAT');


% --------------------------------------------------------------------
function Filtre_adaptatif_Callback(hObject, eventdata, handles)
% hObject    handle to Filtre_adaptatif (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Paramètres du filtre adaptatif
img = handles.courant_data;

% Convertir l'image en niveaux de gris
%grayImg = rgb2gray(img);
grayImg = img;
% Ajouter un bruit gaussien à l'image
noisyImg = imnoise(grayImg, 'gaussian', 0, 0.01);

% Définir la taille de la fenêtre de filtrage
filterSize = 5;
% Définir la variance de la distribution gaussienne pour le calcul du filtre
sigma = 1;

% Appliquer le filtre gaussien à l'image bruitée
gaussFilter = fspecial('gaussian', [filterSize filterSize], sigma);
localMean = imfilter(noisyImg, gaussFilter, 'symmetric');
localVariance = imfilter(noisyImg.^2, gaussFilter, 'symmetric') - localMean.^2;
filteredImg = localMean + (localVariance ./ (localVariance + 0.01)) .* (noisyImg - localMean);

% Afficher les images
imshow(filteredImg);
title('Image filtrée Adaptatif');
