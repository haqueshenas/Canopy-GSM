% Canopy GSM (version 2.1)

%   A code for:
%   Green-gradient based canopy segmentation: A multipurpose image mining model with potential use in crop phenotyping and canopy studies
 
  
%    Authors: Abbas Haghshenas (a), Yahya Emam (b), Saeid Jafarizadeh (C)

%    (a) Shiraz, Iran. 

%    (b) Department of Plant Production and Genetics, Shiraz University, Shiraz, Iran.

%    (c) Former MSc student, Vision Lab, Electrical and Computer Engineering School, Shiraz University, Shiraz, Iran.
    
%   Email: 
%   haqueshenas@gmail.com
%   abbas.haghshenas@shirazu.ac.ir   
           
  
%%  The code is used in the paper entitled:
 
%   Green-gradient based canopy segmentation: A multipurpose image mining model with potential use in crop phenotyping and canopy studies
  
%   (Haghshenas, A., & Emam, Y. (2020). Green-gradient based canopy segmentation: A multipurpose image mining model with potential use in crop
%    phenotyping and canopy studies. Computers and Electronics in Agriculture, 178, 105740. doi:https://doi.org/10.1016/j.compag.2020.105740)
%---------------------------------------------------------------------------------------------------------------

% MIT License

% Copyright (c) 2019-2024 Abbas Haghshenas

% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:

%The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.

% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

%---------------------------------------------------------------------------------------------------------------

% Canopy GSM (version 2.1)

clc
clear
close all

% Setting ST3 variables

% (ST3_P and ST3_m may be set to custom values, depended on the purposes
% and conditions; e.g. in our experiment, they were set to 50 and 20, respectively.)

St3_P = 50;
St3_m = 12;

% Making output subfolders

newSubFolder = sprintf('Graphs', '.\code');
if ~exist(newSubFolder, 'dir')
  mkdir(newSubFolder);
end

newSubFolder = sprintf('Processed images', '.\code');
if ~exist(newSubFolder, 'dir')
  mkdir(newSubFolder);
end

newSubFolder = sprintf('Results', '.\code');
if ~exist(newSubFolder, 'dir')
  mkdir(newSubFolder);
end

% Loading images

disp(['Loading Images from directory (Images)... ']); 
imgpath = './Images/';
imgtype='*.jpg';
images = dir([imgpath imgtype]);
imgtype='*.jpeg';
images = [images ; dir([imgpath imgtype]);]
imgtype='*.tiff';
images = [images ; dir([imgpath imgtype]);]
imgtype='*.tif';
images = [images ; dir([imgpath imgtype]);]

% Table Total_ST2
numChannels = 255;
prefixes = {'R G', 'B G'};
varNames = {'Image'};

for prefixIdx = 1:numel(prefixes)
    prefix = prefixes{prefixIdx};
    
    for channelIdx = 1:numChannels
        varNames{end+1} = [prefix num2str(channelIdx)];
    end
end

Data_total_st2 = cellstr(varNames);

% ---

Data_total_Exp = {'Image','R2_exp_Red','RMSE_exp_Red','exp_a_Red','exp_b_Red','R2_exp_Blue','RMSE_exp_Blue','Exp_a_Blue','Exp_b_Blue'};

Data_total_Num = {'Image'};


for i=2:length(images)+1
Data = {'Green level','Mean Red','Green','Mean Blue',...
'Var Red','Var Blue','Number of pixels','Red_AUC','Blue_AUC','RB_ABC',...
'Red Slope','Blue slope','Red ST3 class','Blue ST3 class','R^2 exp Red','RMSE exp Red ','exp a Red','exp b Red',...
    'R^2 exp Blue','RMSE exp Blue ','exp a Blue','exp b Blue'};

Data4 = {'Image name'};

Image = imread(['./Images/',images(i-1).name]);

R =  (Image(:,:,1));
G =  (Image(:,:,2));
B =  (Image(:,:,3));

disp('=================== New image =======================')
disp(['Processing image(',num2str(i-1),':',images(i-1).name,')... ']);

ImageForgraoundRGB = Image;

BinaryImageRGB = ST1(R,G,B);

ImageForgraound =[];
ImageForgraound(:,:,1) =uint8(ST1(R,G,B)).*Image(:,:,1);
ImageForgraound(:,:,2) =uint8(ST1(R,G,B)).*Image(:,:,2);
ImageForgraound(:,:,3) =uint8(ST1(R,G,B)).*Image(:,:,3);

% Setting background color

r = 150;
g = 0;
b = 150;

ImageForgraound2 = uint8(ChangeBackColor(ImageForgraound,r,g,b));

imwrite(ImageForgraound2,['./Processed images/Processed_',[images(i-1).name],'.jpg'])
 
R1 = R.*uint8(ST1(R,G,B));
G1 = G.*uint8(ST1(R,G,B));
B1 = B.*uint8(ST1(R,G,B));

% GSM (ST2)

for j=1:255
    disp(['Computing level ',num2str(j),' of ',images(i-1).name,'... ']);

Data(j+2,1) = num2cell(j);
Data(j+2,2) = num2cell(mean(mean(R(find(G1==j)))));
Data(j+2,3) = num2cell(j);
Data(j+2,4) = num2cell(mean(mean(B(find(G1==j)))));

%% variance 
R = double(R);
G = double(G);
B = double(B);

Data(j+2,5) = num2cell((var(R(find(G1==j )))));
Data(j+2,6) = num2cell((var(B(find(G1==j )))));
Data(j+2,7) = num2cell(length(R(find(G1==j ))));

tempp = cell2mat(Data(3:end,2));
tempp(isnan(tempp))= 0;
tempp = sum(tempp);

Data(j+2,8) = num2cell(tempp);
tempp = cell2mat(Data(3:end,4));
tempp(isnan(tempp))= 0;
tempp = sum(tempp);

Data(j+2,9) = num2cell(tempp);
Data(j+2,10) = num2cell(cell2mat(Data(j+2,8))-cell2mat(Data(j+2,9)));

end

for j=1:257

    % Curve derivation
    
    if(~isnan(cell2mat(Data(j,2))) & j>2)  
       if(j<257)
            beforind = 1;
            afterind = 1;
            if ~isnan(cell2mat(Data(j,2)))
                temp = cell2mat(Data(j,2));
                while isnan(cell2mat(Data(j+afterind,2))) & j<(j+afterind)-1 
                    afterind=afterind+1;
                end
                while isnan(cell2mat(Data(j-beforind,2))) & j>(j-beforind)+1
                    beforind=beforind+1;
                end
                if (j-beforind)>2 & ~isnan(cell2mat(Data(j-beforind,2)))
                    mat1 = cell2mat(Data(j-beforind,2));
                else
                    beforind = -1;
                    mat1 = cell2mat(Data(j,2));
                end
                if (j+afterind)<=length(Data(:,2)) & ~isnan(cell2mat(Data(j+afterind,2)))
                    mat2 = cell2mat(Data(j+afterind,2));
                else
                    afterind = -1;
                    mat2 = cell2mat(Data(j,2));
                end
                if(afterind==-1)
                    Data(j,11) = num2cell((temp-mat1)/beforind);
                elseif(beforind==-1)
                    Data(j,11) = num2cell((mat2-temp)/afterind);
                else
                    Data(j,11) = num2cell(mean([(mat2-temp)/afterind,(temp-mat1)/beforind]));
                end
                if (j-beforind)>2 & ~isnan(cell2mat(Data(j-beforind,4)))
                    mat1 = cell2mat(Data(j-beforind,4));
                else
                    beforind = -1;
                    mat1 = cell2mat(Data(j,4));
                end
                if (j+afterind)<=length(Data(:,2)) & ~isnan(cell2mat(Data(j+afterind,4)))
                    mat2 = cell2mat(Data(j+afterind,4));
                else
                    afterind = -1;
                    mat2 = cell2mat(Data(j,4));
                end
                temp = cell2mat(Data(j,4));
                if(afterind==-1) 
                    Data(j,12) = num2cell((temp-mat1)/beforind);
                elseif(beforind==-1)
                    Data(j,12) = num2cell((mat2-temp)/afterind);
                else
                    Data(j,12) = num2cell(mean( [(mat2-temp)/afterind,(temp-mat1)/beforind]));
                end
            end
       else
            beforind = 1;
            if ~isnan(cell2mat(Data(j,2)))
                while isnan(cell2mat(Data(j-beforind,2)))
                    beforind=beforind+1;
                end
                if (j-beforind)>2 & ~isnan(cell2mat(Data(j-beforind,2)))
                    mat1 = cell2mat(Data(j-beforind,2));
                else
                     mat1 = cell2mat(Data(j,2));
                end
                mat2 = cell2mat(Data(j,2));
                temp = mat2 - mat1;
                Data(j,11) = num2cell((temp)/(beforind));
                if (j-beforind)>2 & ~isnan(cell2mat(Data(j-beforind,4)))
                    mat1 = cell2mat(Data(j-beforind,4));
                else
                    mat1 = cell2mat(Data(j,4));
                end
                mat2 = cell2mat(Data(j,4));
                temp = mat2 - mat1;
                Data(j,12) = num2cell((temp)/(beforind));
            end
       end
    end
end

% ST3 

temp1 = (ST3(Data(2:end,11),St3_m,St3_P));
temp2 = zeros(1,length(Data(2:end,11))-length(temp1));
temp = [temp2,temp1]';
Data(2:end,13) = num2cell(temp);
temp1 = (ST3(Data(2:end,12),St3_m,St3_P));
temp2 = zeros(1,length(Data(2:end,12))-length(temp1));
temp = [temp2,temp1]';
Data(2:end,14) = num2cell(temp);

% Fitting exponential equations

[R1,R2,X_,Y_] = ExpFitting(Data(3:end,1),Data(3:end,2));
 Data(2,15) = num2cell(R2.rsquare);
 Data(2,16) = num2cell(R2.rmse);
 Data(2,17) = num2cell(R1.a);
 Data(2,18) = num2cell(R1.b);
[R1,R2,X_,Y_] = ExpFitting(Data(3:end,3),Data(3:end,4));
 Data(2,19) = num2cell(R2.rsquare);
 Data(2,20) = num2cell(R2.rmse);
 Data(2,21) = num2cell(R1.a);
 Data(2,22) = num2cell(R1.b); 

% Writing single tables

T = table(Data(2:end,1), Data(2:end,2), Data(2:end,3), Data(2:end,4), Data(2:end,5),...
    Data(2:end,6), Data(2:end,7),Data(2:end,8), Data(2:end,9), Data(2:end,10),Data(2:end,11),Data(2:end,12),Data(2:end,13),Data(2:end,14),...
    Data(2:end,15),Data(2:end,16),Data(2:end,17),Data(2:end,18),Data(2:end,19),Data(2:end,20),Data(2:end,21),Data(2:end,22));

T.Properties.VariableNames={'Green_level','Mean_Red','Green','Mean_Blue','Var_Red','Var_Blue','Number_of_pixels','Red_AUC','Blue_AUC','RB_ABC'...
    ,'Red_Slope','Blue_slope','Red_ST3_class','Blue_ST3_class','R2_exp_Red','RMSE_exp_Red ','exp_a_Red','exp_b_Red',...
    'R2_exp_Blue','RMSE_exp_Blue ','exp_a_Blue','exp_b_Blue'};

writetable(T,['./Results/Results of image ',[images(i-1).name],'.csv'],'Delimiter',',');

Data_total_st2(i,1) = cellstr(images(i-1).name);

for ii=2:length(Data(:,2))-1
    Data_total_st2(i,ii)=Data(ii+1,2);
end

for ii=2:length(Data(:,4))-1
    Data_total_st2(i,ii+255)=Data(ii+1,4);
end


Data_total_num(i,1) = cellstr(images(i-1).name);

for ii=2:length(Data(:,7))-1
    Data_total_num(i,ii)=Data(ii+1,7);
end

Data_total_Exp(i,1) = cellstr(images(i-1).name);
Data_total_Exp(i,2) = Data(2,15);
Data_total_Exp(i,3) = Data(2,16);
Data_total_Exp(i,4) = Data(2,17);
Data_total_Exp(i,5) = Data(2,18);
Data_total_Exp(i,6) = Data(2,19);
Data_total_Exp(i,7) = Data(2,20);
Data_total_Exp(i,8) = Data(2,21);
Data_total_Exp(i,9) = Data(2,22);

GSMRow = zeros(length(Data(3:end,1)),1);

for k=3:length(Data(:,1))
    GSMRow(k,1) = cell2mat(Data(k,1)); 
end

GSMMean_Red = zeros(length(Data(:,2)),1);

for k=3:length(Data(:,2))
    GSMMean_Red(k,1) = cell2mat(Data(k,2)); 
end

GSMGreen = zeros(length(Data(:,3)),1);

for k=3:length(Data(:,3))
    GSMGreen(k,1) = cell2mat(Data(k,3)); 
end

GSMMean_Blue = zeros(length(Data(:,4)),1);

for k=3:length(Data(:,4))
    GSMMean_Blue(k,1) = cell2mat(Data(k,4)); 
end

GSMNO_of_Pixel = [];

for k=3:length(Data(:,7))
    GSMNO_of_Pixel = [GSMNO_of_Pixel;cell2mat(Data(k,7))]; 
end


% Graphs

h = figure();
ax1 = subplot(2,1,1);
hold on
xlim manual
xlim([-1 256])
ylim manual
ylim([-1 256])
plot(GSMRow,GSMMean_Red,'.r');
plot(GSMRow,GSMGreen,'.g');
plot(GSMRow,GSMMean_Blue,'.b');
hold off
title(ax1,'GSM Graph');
ylabel(ax1,'Mean Red, Green, Mean Blue (0-255)');
xlabel(ax1,'Green Level (1-255)');
pbaspect([1 1 1])
ax2 = subplot(2,1,2);
hold on
xlim manual
xlim([-1 256])
ylim auto
plot(GSMRow(3:end,1),GSMNO_of_Pixel,'.k');
title(ax2,'Number of pixels');
ylabel(ax2,'Number of pixels');
xlabel(ax2,'Green Level (1-255)');
pbaspect([1 1 1])
hold off
set(gcf, 'PaperUnits', 'inches');
x_width=7.5; y_width=9.25;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);
saveas(h,['./Graphs/Graphs of image ',[images(i-1).name],'.pdf'],'pdf')
close all
end

T_Exp = table(Data_total_Exp(2:end,1),Data_total_Exp(2:end,2),Data_total_Exp(2:end,3),Data_total_Exp(2:end,4),Data_total_Exp(2:end,5),...
Data_total_Exp(2:end,6),Data_total_Exp(2:end,7),Data_total_Exp(2:end,8),Data_total_Exp(2:end,9));

T_Exp.Properties.VariableNames={'Image','R2_exp_Red','RMSE_exp_Red','Exp_a_Red','Exp_b_Red','R2_exp_Blue','RMSE_exp_Blue','Exp_a_Blue','Exp_b_Blue'};
writetable(T_Exp,'./Results/Total GSM Exp.csv','Delimiter',',');


T_St2 = table(Data_total_st2( 2:end ,1),Data_total_st2( 2:end ,2),Data_total_st2( 2:end ,3),Data_total_st2( 2:end ,4),Data_total_st2( 2:end ,5),Data_total_st2( 2:end ,6),Data_total_st2( 2:end ,7),...
Data_total_st2( 2:end ,8),Data_total_st2( 2:end ,9),Data_total_st2( 2:end ,10),Data_total_st2( 2:end ,11),Data_total_st2( 2:end ,12),Data_total_st2( 2:end ,13),...
Data_total_st2( 2:end ,14),Data_total_st2( 2:end ,15),Data_total_st2( 2:end ,16),Data_total_st2( 2:end ,17),Data_total_st2( 2:end ,18),Data_total_st2( 2:end ,19),...
Data_total_st2( 2:end ,20),Data_total_st2( 2:end ,21),Data_total_st2( 2:end ,22),Data_total_st2( 2:end ,23),Data_total_st2( 2:end ,24),Data_total_st2( 2:end ,25),...
Data_total_st2( 2:end ,26),Data_total_st2( 2:end ,27),Data_total_st2( 2:end ,28),Data_total_st2( 2:end ,29),Data_total_st2( 2:end ,30),Data_total_st2( 2:end ,31),...
Data_total_st2( 2:end ,32),Data_total_st2( 2:end ,33),Data_total_st2( 2:end ,34),Data_total_st2( 2:end ,35),Data_total_st2( 2:end ,36),Data_total_st2( 2:end ,37),...
Data_total_st2( 2:end ,38),Data_total_st2( 2:end ,39),Data_total_st2( 2:end ,40),Data_total_st2( 2:end ,41),Data_total_st2( 2:end ,42),Data_total_st2( 2:end ,43),...
Data_total_st2( 2:end ,44),Data_total_st2( 2:end ,45),Data_total_st2( 2:end ,46),Data_total_st2( 2:end ,47),Data_total_st2( 2:end ,48),Data_total_st2( 2:end ,49),...
Data_total_st2( 2:end ,50),Data_total_st2( 2:end ,51),Data_total_st2( 2:end ,52),Data_total_st2( 2:end ,53),Data_total_st2( 2:end ,54),Data_total_st2( 2:end ,55),...
Data_total_st2( 2:end ,56),Data_total_st2( 2:end ,57),Data_total_st2( 2:end ,58),Data_total_st2( 2:end ,59),Data_total_st2( 2:end ,60),Data_total_st2( 2:end ,61),...
Data_total_st2( 2:end ,62),Data_total_st2( 2:end ,63),Data_total_st2( 2:end ,64),Data_total_st2( 2:end ,65),Data_total_st2( 2:end ,66),Data_total_st2( 2:end ,67),...
Data_total_st2( 2:end ,68),Data_total_st2( 2:end ,69),Data_total_st2( 2:end ,70),Data_total_st2( 2:end ,71),Data_total_st2( 2:end ,72),Data_total_st2( 2:end ,73),...
Data_total_st2( 2:end ,74),Data_total_st2( 2:end ,75),Data_total_st2( 2:end ,76),Data_total_st2( 2:end ,77),Data_total_st2( 2:end ,78),Data_total_st2( 2:end ,79),...
Data_total_st2( 2:end ,80),Data_total_st2( 2:end ,81),Data_total_st2( 2:end ,82),Data_total_st2( 2:end ,83),Data_total_st2( 2:end ,84),Data_total_st2( 2:end ,85),...
Data_total_st2( 2:end ,86),Data_total_st2( 2:end ,87),Data_total_st2( 2:end ,88),Data_total_st2( 2:end ,89),Data_total_st2( 2:end ,90),Data_total_st2( 2:end ,91),...
Data_total_st2( 2:end ,92),Data_total_st2( 2:end ,93),Data_total_st2( 2:end ,94),Data_total_st2( 2:end ,95),Data_total_st2( 2:end ,96),Data_total_st2( 2:end ,97),...
Data_total_st2( 2:end ,98),Data_total_st2( 2:end ,99),Data_total_st2( 2:end ,100),Data_total_st2( 2:end ,101),Data_total_st2( 2:end ,102),Data_total_st2( 2:end ,103),...
Data_total_st2( 2:end ,104),Data_total_st2( 2:end ,105),Data_total_st2( 2:end ,106),Data_total_st2( 2:end ,107),Data_total_st2( 2:end ,108),Data_total_st2( 2:end ,109),...
Data_total_st2( 2:end ,110),Data_total_st2( 2:end ,111),Data_total_st2( 2:end ,112),Data_total_st2( 2:end ,113),Data_total_st2( 2:end ,114),Data_total_st2( 2:end ,115),...
Data_total_st2( 2:end ,116),Data_total_st2( 2:end ,117),Data_total_st2( 2:end ,118),Data_total_st2( 2:end ,119),Data_total_st2( 2:end ,120),Data_total_st2( 2:end ,121),...
Data_total_st2( 2:end ,122),Data_total_st2( 2:end ,123),Data_total_st2( 2:end ,124),Data_total_st2( 2:end ,125),Data_total_st2( 2:end ,126),Data_total_st2( 2:end ,127),...
Data_total_st2( 2:end ,128),Data_total_st2( 2:end ,129),Data_total_st2( 2:end ,130),Data_total_st2( 2:end ,131),Data_total_st2( 2:end ,132),Data_total_st2( 2:end ,133),...
Data_total_st2( 2:end ,134),Data_total_st2( 2:end ,135),Data_total_st2( 2:end ,136),Data_total_st2( 2:end ,137),Data_total_st2( 2:end ,138),Data_total_st2( 2:end ,139),...
Data_total_st2( 2:end ,140),Data_total_st2( 2:end ,141),Data_total_st2( 2:end ,142),Data_total_st2( 2:end ,143),Data_total_st2( 2:end ,144),Data_total_st2( 2:end ,145),...
Data_total_st2( 2:end ,146),Data_total_st2( 2:end ,147),Data_total_st2( 2:end ,148),Data_total_st2( 2:end ,149),Data_total_st2( 2:end ,150),Data_total_st2( 2:end ,151),...
Data_total_st2( 2:end ,152),Data_total_st2( 2:end ,153),Data_total_st2( 2:end ,154),Data_total_st2( 2:end ,155),Data_total_st2( 2:end ,156),Data_total_st2( 2:end ,157),...
Data_total_st2( 2:end ,158),Data_total_st2( 2:end ,159),Data_total_st2( 2:end ,160),Data_total_st2( 2:end ,161),Data_total_st2( 2:end ,162),Data_total_st2( 2:end ,163),...
Data_total_st2( 2:end ,164),Data_total_st2( 2:end ,165),Data_total_st2( 2:end ,166),Data_total_st2( 2:end ,167),Data_total_st2( 2:end ,168),Data_total_st2( 2:end ,169),...
Data_total_st2( 2:end ,170),Data_total_st2( 2:end ,171),Data_total_st2( 2:end ,172),Data_total_st2( 2:end ,173),Data_total_st2( 2:end ,174),Data_total_st2( 2:end ,175),...
Data_total_st2( 2:end ,176),Data_total_st2( 2:end ,177),Data_total_st2( 2:end ,178),Data_total_st2( 2:end ,179),Data_total_st2( 2:end ,180),Data_total_st2( 2:end ,181),...
Data_total_st2( 2:end ,182),Data_total_st2( 2:end ,183),Data_total_st2( 2:end ,184),Data_total_st2( 2:end ,185),Data_total_st2( 2:end ,186),Data_total_st2( 2:end ,187),...
Data_total_st2( 2:end ,188),Data_total_st2( 2:end ,189),Data_total_st2( 2:end ,190),Data_total_st2( 2:end ,191),Data_total_st2( 2:end ,192),Data_total_st2( 2:end ,193),...
Data_total_st2( 2:end ,194),Data_total_st2( 2:end ,195),Data_total_st2( 2:end ,196),Data_total_st2( 2:end ,197),Data_total_st2( 2:end ,198),Data_total_st2( 2:end ,199),...
Data_total_st2( 2:end ,200),Data_total_st2( 2:end ,201),Data_total_st2( 2:end ,202),Data_total_st2( 2:end ,203),Data_total_st2( 2:end ,204),Data_total_st2( 2:end ,205),...
Data_total_st2( 2:end ,206),Data_total_st2( 2:end ,207),Data_total_st2( 2:end ,208),Data_total_st2( 2:end ,209),Data_total_st2( 2:end ,210),Data_total_st2( 2:end ,211),...
Data_total_st2( 2:end ,212),Data_total_st2( 2:end ,213),Data_total_st2( 2:end ,214),Data_total_st2( 2:end ,215),Data_total_st2( 2:end ,216),Data_total_st2( 2:end ,217),...
Data_total_st2( 2:end ,218),Data_total_st2( 2:end ,219),Data_total_st2( 2:end ,220),Data_total_st2( 2:end ,221),Data_total_st2( 2:end ,222),Data_total_st2( 2:end ,223),...
Data_total_st2( 2:end ,224),Data_total_st2( 2:end ,225),Data_total_st2( 2:end ,226),Data_total_st2( 2:end ,227),Data_total_st2( 2:end ,228),Data_total_st2( 2:end ,229),...
Data_total_st2( 2:end ,230),Data_total_st2( 2:end ,231),Data_total_st2( 2:end ,232),Data_total_st2( 2:end ,233),Data_total_st2( 2:end ,234),Data_total_st2( 2:end ,235),...
Data_total_st2( 2:end ,236),Data_total_st2( 2:end ,237),Data_total_st2( 2:end ,238),Data_total_st2( 2:end ,239),Data_total_st2( 2:end ,240),Data_total_st2( 2:end ,241),...
Data_total_st2( 2:end ,242),Data_total_st2( 2:end ,243),Data_total_st2( 2:end ,244),Data_total_st2( 2:end ,245),Data_total_st2( 2:end ,246),Data_total_st2( 2:end ,247),...
Data_total_st2( 2:end ,248),Data_total_st2( 2:end ,249),Data_total_st2( 2:end ,250),Data_total_st2( 2:end ,251),Data_total_st2( 2:end ,252),Data_total_st2( 2:end ,253),...
Data_total_st2( 2:end ,254),Data_total_st2( 2:end ,255),Data_total_st2( 2:end ,256),Data_total_st2( 2:end ,257),Data_total_st2( 2:end ,258),Data_total_st2( 2:end ,259),...
Data_total_st2( 2:end ,260),Data_total_st2( 2:end ,261),Data_total_st2( 2:end ,262),Data_total_st2( 2:end ,263),Data_total_st2( 2:end ,264),Data_total_st2( 2:end ,265),...
Data_total_st2( 2:end ,266),Data_total_st2( 2:end ,267),Data_total_st2( 2:end ,268),Data_total_st2( 2:end ,269),Data_total_st2( 2:end ,270),Data_total_st2( 2:end ,271),...
Data_total_st2( 2:end ,272),Data_total_st2( 2:end ,273),Data_total_st2( 2:end ,274),Data_total_st2( 2:end ,275),Data_total_st2( 2:end ,276),Data_total_st2( 2:end ,277),...
Data_total_st2( 2:end ,278),Data_total_st2( 2:end ,279),Data_total_st2( 2:end ,280),Data_total_st2( 2:end ,281),Data_total_st2( 2:end ,282),Data_total_st2( 2:end ,283),...
Data_total_st2( 2:end ,284),Data_total_st2( 2:end ,285),Data_total_st2( 2:end ,286),Data_total_st2( 2:end ,287),Data_total_st2( 2:end ,288),Data_total_st2( 2:end ,289),...
Data_total_st2( 2:end ,290),Data_total_st2( 2:end ,291),Data_total_st2( 2:end ,292),Data_total_st2( 2:end ,293),Data_total_st2( 2:end ,294),Data_total_st2( 2:end ,295),...
Data_total_st2( 2:end ,296),Data_total_st2( 2:end ,297),Data_total_st2( 2:end ,298),Data_total_st2( 2:end ,299),Data_total_st2( 2:end ,300),Data_total_st2( 2:end ,301),...
Data_total_st2( 2:end ,302),Data_total_st2( 2:end ,303),Data_total_st2( 2:end ,304),Data_total_st2( 2:end ,305),Data_total_st2( 2:end ,306),Data_total_st2( 2:end ,307),...
Data_total_st2( 2:end ,308),Data_total_st2( 2:end ,309),Data_total_st2( 2:end ,310),Data_total_st2( 2:end ,311),Data_total_st2( 2:end ,312),Data_total_st2( 2:end ,313),...
Data_total_st2( 2:end ,314),Data_total_st2( 2:end ,315),Data_total_st2( 2:end ,316),Data_total_st2( 2:end ,317),Data_total_st2( 2:end ,318),Data_total_st2( 2:end ,319),...
Data_total_st2( 2:end ,320),Data_total_st2( 2:end ,321),Data_total_st2( 2:end ,322),Data_total_st2( 2:end ,323),Data_total_st2( 2:end ,324),Data_total_st2( 2:end ,325),...
Data_total_st2( 2:end ,326),Data_total_st2( 2:end ,327),Data_total_st2( 2:end ,328),Data_total_st2( 2:end ,329),Data_total_st2( 2:end ,330),Data_total_st2( 2:end ,331),...
Data_total_st2( 2:end ,332),Data_total_st2( 2:end ,333),Data_total_st2( 2:end ,334),Data_total_st2( 2:end ,335),Data_total_st2( 2:end ,336),Data_total_st2( 2:end ,337),...
Data_total_st2( 2:end ,338),Data_total_st2( 2:end ,339),Data_total_st2( 2:end ,340),Data_total_st2( 2:end ,341),Data_total_st2( 2:end ,342),Data_total_st2( 2:end ,343),...
Data_total_st2( 2:end ,344),Data_total_st2( 2:end ,345),Data_total_st2( 2:end ,346),Data_total_st2( 2:end ,347),Data_total_st2( 2:end ,348),Data_total_st2( 2:end ,349),...
Data_total_st2( 2:end ,350),Data_total_st2( 2:end ,351),Data_total_st2( 2:end ,352),Data_total_st2( 2:end ,353),Data_total_st2( 2:end ,354),Data_total_st2( 2:end ,355),...
Data_total_st2( 2:end ,356),Data_total_st2( 2:end ,357),Data_total_st2( 2:end ,358),Data_total_st2( 2:end ,359),Data_total_st2( 2:end ,360),Data_total_st2( 2:end ,361),...
Data_total_st2( 2:end ,362),Data_total_st2( 2:end ,363),Data_total_st2( 2:end ,364),Data_total_st2( 2:end ,365),Data_total_st2( 2:end ,366),Data_total_st2( 2:end ,367),...
Data_total_st2( 2:end ,368),Data_total_st2( 2:end ,369),Data_total_st2( 2:end ,370),Data_total_st2( 2:end ,371),Data_total_st2( 2:end ,372),Data_total_st2( 2:end ,373),...
Data_total_st2( 2:end ,374),Data_total_st2( 2:end ,375),Data_total_st2( 2:end ,376),Data_total_st2( 2:end ,377),Data_total_st2( 2:end ,378),Data_total_st2( 2:end ,379),...
Data_total_st2( 2:end ,380),Data_total_st2( 2:end ,381),Data_total_st2( 2:end ,382),Data_total_st2( 2:end ,383),Data_total_st2( 2:end ,384),Data_total_st2( 2:end ,385),...
Data_total_st2( 2:end ,386),Data_total_st2( 2:end ,387),Data_total_st2( 2:end ,388),Data_total_st2( 2:end ,389),Data_total_st2( 2:end ,390),Data_total_st2( 2:end ,391),...
Data_total_st2( 2:end ,392),Data_total_st2( 2:end ,393),Data_total_st2( 2:end ,394),Data_total_st2( 2:end ,395),Data_total_st2( 2:end ,396),Data_total_st2( 2:end ,397),...
Data_total_st2( 2:end ,398),Data_total_st2( 2:end ,399),Data_total_st2( 2:end ,400),Data_total_st2( 2:end ,401),Data_total_st2( 2:end ,402),Data_total_st2( 2:end ,403),...
Data_total_st2( 2:end ,404),Data_total_st2( 2:end ,405),Data_total_st2( 2:end ,406),Data_total_st2( 2:end ,407),Data_total_st2( 2:end ,408),Data_total_st2( 2:end ,409),...
Data_total_st2( 2:end ,410),Data_total_st2( 2:end ,411),Data_total_st2( 2:end ,412),Data_total_st2( 2:end ,413),Data_total_st2( 2:end ,414),Data_total_st2( 2:end ,415),...
Data_total_st2( 2:end ,416),Data_total_st2( 2:end ,417),Data_total_st2( 2:end ,418),Data_total_st2( 2:end ,419),Data_total_st2( 2:end ,420),Data_total_st2( 2:end ,421),...
Data_total_st2( 2:end ,422),Data_total_st2( 2:end ,423),Data_total_st2( 2:end ,424),Data_total_st2( 2:end ,425),Data_total_st2( 2:end ,426),Data_total_st2( 2:end ,427),...
Data_total_st2( 2:end ,428),Data_total_st2( 2:end ,429),Data_total_st2( 2:end ,430),Data_total_st2( 2:end ,431),Data_total_st2( 2:end ,432),Data_total_st2( 2:end ,433),...
Data_total_st2( 2:end ,434),Data_total_st2( 2:end ,435),Data_total_st2( 2:end ,436),Data_total_st2( 2:end ,437),Data_total_st2( 2:end ,438),Data_total_st2( 2:end ,439),...
Data_total_st2( 2:end ,440),Data_total_st2( 2:end ,441),Data_total_st2( 2:end ,442),Data_total_st2( 2:end ,443),Data_total_st2( 2:end ,444),Data_total_st2( 2:end ,445),...
Data_total_st2( 2:end ,446),Data_total_st2( 2:end ,447),Data_total_st2( 2:end ,448),Data_total_st2( 2:end ,449),Data_total_st2( 2:end ,450),Data_total_st2( 2:end ,451),...
Data_total_st2( 2:end ,452),Data_total_st2( 2:end ,453),Data_total_st2( 2:end ,454),Data_total_st2( 2:end ,455),Data_total_st2( 2:end ,456),Data_total_st2( 2:end ,457),...
Data_total_st2( 2:end ,458),Data_total_st2( 2:end ,459),Data_total_st2( 2:end ,460),Data_total_st2( 2:end ,461),Data_total_st2( 2:end ,462),Data_total_st2( 2:end ,463),...
Data_total_st2( 2:end ,464),Data_total_st2( 2:end ,465),Data_total_st2( 2:end ,466),Data_total_st2( 2:end ,467),Data_total_st2( 2:end ,468),Data_total_st2( 2:end ,469),...
Data_total_st2( 2:end ,470),Data_total_st2( 2:end ,471),Data_total_st2( 2:end ,472),Data_total_st2( 2:end ,473),Data_total_st2( 2:end ,474),Data_total_st2( 2:end ,475),...
Data_total_st2( 2:end ,476),Data_total_st2( 2:end ,477),Data_total_st2( 2:end ,478),Data_total_st2( 2:end ,479),Data_total_st2( 2:end ,480),Data_total_st2( 2:end ,481),...
Data_total_st2( 2:end ,482),Data_total_st2( 2:end ,483),Data_total_st2( 2:end ,484),Data_total_st2( 2:end ,485),Data_total_st2( 2:end ,486),Data_total_st2( 2:end ,487),...
Data_total_st2( 2:end ,488),Data_total_st2( 2:end ,489),Data_total_st2( 2:end ,490),Data_total_st2( 2:end ,491),Data_total_st2( 2:end ,492),Data_total_st2( 2:end ,493),...
Data_total_st2( 2:end ,494),Data_total_st2( 2:end ,495),Data_total_st2( 2:end ,496),Data_total_st2( 2:end ,497),Data_total_st2( 2:end ,498),Data_total_st2( 2:end ,499),...
Data_total_st2( 2:end ,500),Data_total_st2( 2:end ,501),Data_total_st2( 2:end ,502),Data_total_st2( 2:end ,503),Data_total_st2( 2:end ,504),Data_total_st2( 2:end ,505),...
Data_total_st2( 2:end ,506),Data_total_st2( 2:end ,507),Data_total_st2( 2:end ,508),Data_total_st2( 2:end ,509),Data_total_st2( 2:end ,510),Data_total_st2( 2:end ,511));

%Table ST2 Headers

numChannels = 255;
prefixes = {'R_G', 'B_G'};
varNames = {'Image'};

for prefixIdx = 1:numel(prefixes)
    prefix = prefixes{prefixIdx};
    
    for channelIdx = 1:numChannels
        varNames{end+1} = [prefix num2str(channelIdx)];
    end
end

T_St2.Properties.VariableNames = varNames;

% -----

writetable(T_St2,'./Results/Total GSM ST2.csv','Delimiter',',');

T_Num = table(Data_total_num( 2:end ,1),Data_total_num( 2:end ,2),Data_total_num( 2:end ,3),Data_total_num( 2:end ,4),Data_total_num( 2:end ,5),Data_total_num( 2:end ,6),Data_total_num( 2:end ,7),...
Data_total_num( 2:end ,8),Data_total_num( 2:end ,9),Data_total_num( 2:end ,10),Data_total_num( 2:end ,11),Data_total_num( 2:end ,12),Data_total_num( 2:end ,13),...
Data_total_num( 2:end ,14),Data_total_num( 2:end ,15),Data_total_num( 2:end ,16),Data_total_num( 2:end ,17),Data_total_num( 2:end ,18),Data_total_num( 2:end ,19),...
Data_total_num( 2:end ,20),Data_total_num( 2:end ,21),Data_total_num( 2:end ,22),Data_total_num( 2:end ,23),Data_total_num( 2:end ,24),Data_total_num( 2:end ,25),...
Data_total_num( 2:end ,26),Data_total_num( 2:end ,27),Data_total_num( 2:end ,28),Data_total_num( 2:end ,29),Data_total_num( 2:end ,30),Data_total_num( 2:end ,31),...
Data_total_num( 2:end ,32),Data_total_num( 2:end ,33),Data_total_num( 2:end ,34),Data_total_num( 2:end ,35),Data_total_num( 2:end ,36),Data_total_num( 2:end ,37),...
Data_total_num( 2:end ,38),Data_total_num( 2:end ,39),Data_total_num( 2:end ,40),Data_total_num( 2:end ,41),Data_total_num( 2:end ,42),Data_total_num( 2:end ,43),...
Data_total_num( 2:end ,44),Data_total_num( 2:end ,45),Data_total_num( 2:end ,46),Data_total_num( 2:end ,47),Data_total_num( 2:end ,48),Data_total_num( 2:end ,49),...
Data_total_num( 2:end ,50),Data_total_num( 2:end ,51),Data_total_num( 2:end ,52),Data_total_num( 2:end ,53),Data_total_num( 2:end ,54),Data_total_num( 2:end ,55),...
Data_total_num( 2:end ,56),Data_total_num( 2:end ,57),Data_total_num( 2:end ,58),Data_total_num( 2:end ,59),Data_total_num( 2:end ,60),Data_total_num( 2:end ,61),...
Data_total_num( 2:end ,62),Data_total_num( 2:end ,63),Data_total_num( 2:end ,64),Data_total_num( 2:end ,65),Data_total_num( 2:end ,66),Data_total_num( 2:end ,67),...
Data_total_num( 2:end ,68),Data_total_num( 2:end ,69),Data_total_num( 2:end ,70),Data_total_num( 2:end ,71),Data_total_num( 2:end ,72),Data_total_num( 2:end ,73),...
Data_total_num( 2:end ,74),Data_total_num( 2:end ,75),Data_total_num( 2:end ,76),Data_total_num( 2:end ,77),Data_total_num( 2:end ,78),Data_total_num( 2:end ,79),...
Data_total_num( 2:end ,80),Data_total_num( 2:end ,81),Data_total_num( 2:end ,82),Data_total_num( 2:end ,83),Data_total_num( 2:end ,84),Data_total_num( 2:end ,85),...
Data_total_num( 2:end ,86),Data_total_num( 2:end ,87),Data_total_num( 2:end ,88),Data_total_num( 2:end ,89),Data_total_num( 2:end ,90),Data_total_num( 2:end ,91),...
Data_total_num( 2:end ,92),Data_total_num( 2:end ,93),Data_total_num( 2:end ,94),Data_total_num( 2:end ,95),Data_total_num( 2:end ,96),Data_total_num( 2:end ,97),...
Data_total_num( 2:end ,98),Data_total_num( 2:end ,99),Data_total_num( 2:end ,100),Data_total_num( 2:end ,101),Data_total_num( 2:end ,102),Data_total_num( 2:end ,103),...
Data_total_num( 2:end ,104),Data_total_num( 2:end ,105),Data_total_num( 2:end ,106),Data_total_num( 2:end ,107),Data_total_num( 2:end ,108),Data_total_num( 2:end ,109),...
Data_total_num( 2:end ,110),Data_total_num( 2:end ,111),Data_total_num( 2:end ,112),Data_total_num( 2:end ,113),Data_total_num( 2:end ,114),Data_total_num( 2:end ,115),...
Data_total_num( 2:end ,116),Data_total_num( 2:end ,117),Data_total_num( 2:end ,118),Data_total_num( 2:end ,119),Data_total_num( 2:end ,120),Data_total_num( 2:end ,121),...
Data_total_num( 2:end ,122),Data_total_num( 2:end ,123),Data_total_num( 2:end ,124),Data_total_num( 2:end ,125),Data_total_num( 2:end ,126),Data_total_num( 2:end ,127),...
Data_total_num( 2:end ,128),Data_total_num( 2:end ,129),Data_total_num( 2:end ,130),Data_total_num( 2:end ,131),Data_total_num( 2:end ,132),Data_total_num( 2:end ,133),...
Data_total_num( 2:end ,134),Data_total_num( 2:end ,135),Data_total_num( 2:end ,136),Data_total_num( 2:end ,137),Data_total_num( 2:end ,138),Data_total_num( 2:end ,139),...
Data_total_num( 2:end ,140),Data_total_num( 2:end ,141),Data_total_num( 2:end ,142),Data_total_num( 2:end ,143),Data_total_num( 2:end ,144),Data_total_num( 2:end ,145),...
Data_total_num( 2:end ,146),Data_total_num( 2:end ,147),Data_total_num( 2:end ,148),Data_total_num( 2:end ,149),Data_total_num( 2:end ,150),Data_total_num( 2:end ,151),...
Data_total_num( 2:end ,152),Data_total_num( 2:end ,153),Data_total_num( 2:end ,154),Data_total_num( 2:end ,155),Data_total_num( 2:end ,156),Data_total_num( 2:end ,157),...
Data_total_num( 2:end ,158),Data_total_num( 2:end ,159),Data_total_num( 2:end ,160),Data_total_num( 2:end ,161),Data_total_num( 2:end ,162),Data_total_num( 2:end ,163),...
Data_total_num( 2:end ,164),Data_total_num( 2:end ,165),Data_total_num( 2:end ,166),Data_total_num( 2:end ,167),Data_total_num( 2:end ,168),Data_total_num( 2:end ,169),...
Data_total_num( 2:end ,170),Data_total_num( 2:end ,171),Data_total_num( 2:end ,172),Data_total_num( 2:end ,173),Data_total_num( 2:end ,174),Data_total_num( 2:end ,175),...
Data_total_num( 2:end ,176),Data_total_num( 2:end ,177),Data_total_num( 2:end ,178),Data_total_num( 2:end ,179),Data_total_num( 2:end ,180),Data_total_num( 2:end ,181),...
Data_total_num( 2:end ,182),Data_total_num( 2:end ,183),Data_total_num( 2:end ,184),Data_total_num( 2:end ,185),Data_total_num( 2:end ,186),Data_total_num( 2:end ,187),...
Data_total_num( 2:end ,188),Data_total_num( 2:end ,189),Data_total_num( 2:end ,190),Data_total_num( 2:end ,191),Data_total_num( 2:end ,192),Data_total_num( 2:end ,193),...
Data_total_num( 2:end ,194),Data_total_num( 2:end ,195),Data_total_num( 2:end ,196),Data_total_num( 2:end ,197),Data_total_num( 2:end ,198),Data_total_num( 2:end ,199),...
Data_total_num( 2:end ,200),Data_total_num( 2:end ,201),Data_total_num( 2:end ,202),Data_total_num( 2:end ,203),Data_total_num( 2:end ,204),Data_total_num( 2:end ,205),...
Data_total_num( 2:end ,206),Data_total_num( 2:end ,207),Data_total_num( 2:end ,208),Data_total_num( 2:end ,209),Data_total_num( 2:end ,210),Data_total_num( 2:end ,211),...
Data_total_num( 2:end ,212),Data_total_num( 2:end ,213),Data_total_num( 2:end ,214),Data_total_num( 2:end ,215),Data_total_num( 2:end ,216),Data_total_num( 2:end ,217),...
Data_total_num( 2:end ,218),Data_total_num( 2:end ,219),Data_total_num( 2:end ,220),Data_total_num( 2:end ,221),Data_total_num( 2:end ,222),Data_total_num( 2:end ,223),...
Data_total_num( 2:end ,224),Data_total_num( 2:end ,225),Data_total_num( 2:end ,226),Data_total_num( 2:end ,227),Data_total_num( 2:end ,228),Data_total_num( 2:end ,229),...
Data_total_num( 2:end ,230),Data_total_num( 2:end ,231),Data_total_num( 2:end ,232),Data_total_num( 2:end ,233),Data_total_num( 2:end ,234),Data_total_num( 2:end ,235),...
Data_total_num( 2:end ,236),Data_total_num( 2:end ,237),Data_total_num( 2:end ,238),Data_total_num( 2:end ,239),Data_total_num( 2:end ,240),Data_total_num( 2:end ,241),...
Data_total_num( 2:end ,242),Data_total_num( 2:end ,243),Data_total_num( 2:end ,244),Data_total_num( 2:end ,245),Data_total_num( 2:end ,246),Data_total_num( 2:end ,247),...
Data_total_num( 2:end ,248),Data_total_num( 2:end ,249),Data_total_num( 2:end ,250),Data_total_num( 2:end ,251),Data_total_num( 2:end ,252),Data_total_num( 2:end ,253),...
Data_total_num( 2:end ,254),Data_total_num( 2:end ,255), Data_total_num( 2:end ,256));

% Table T_Num Headers

T_Num.Properties.VariableNames = [{'Image'}, strcat('G', arrayfun(@num2str, 1:255, 'UniformOutput', false))];


writetable(T_Num,'./Results/Total GSM Num.csv','Delimiter',',');


disp('===== Processing completed! =====')

% End