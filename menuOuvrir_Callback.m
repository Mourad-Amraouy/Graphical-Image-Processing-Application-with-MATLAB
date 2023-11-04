
% --------------------------------------------------------------------
function menuOuvrir_Callback(hObject, eventdata, handles)
% hObject    handle to menuOuvrir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file,path] = uigetfile('*.*');
handles.ima = imread(sprintf('%s',path,file));
axes(handles.imgO)
handles.courant_data = handles.ima;
subimage(handles.courant_data);
axes(handles.imgT)
% handles.ima_traite = 0;
% subimage(handles.ima_traite);

subimage(handles.courant_data);

handles.output = hObject;
guidata(hObject, handles);
% --------------------------------------------------------------------
function menuEnregistrer_Callback(hObject, eventdata, handles)
% hObject    handle to menuEnregistrer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
image = handles.courant_data;
[file,path] = uiputfile('*.png','Enregistrer Votre Image ...');
imwrite(image, sprintf('%s',path,file),'png');

% --------------------------------------------------------------------
function menuQuitter_Callback(hObject, eventdata, handles)
% hObject    handle to menuQuitter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1)

% --------------------------------------------------------------------
function menuMedian_Callback(hObject, eventdata, handles)
% hObject    handle to menuMedian (see GCBO)
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
