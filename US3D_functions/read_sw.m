function [PXDATA,RES,IM,LANDMARK,OBJECT] = read_sw(sw_filename)
%%READ_SW Reads some parameters from Stradwin sw files.
if nargin == 0
    [file,path] = uigetfile('*.sw');
    if file == 0;return;end
    sw_filename = fullfile(path,file);
end
fid = fopen(sw_filename);
tline = fgetl(fid);
r=0;i=0;l=0;o=0;
LANDMARK=[];
OBJECT = [];
IM = [];
RED=[];
while ischar(tline)
    if startsWith(tline,'IM') == true
        i = i + 1;
        str = textscan(tline,'%s %f %f %f %f %f %f %f');
%         IM(i,1:6) = cell2mat(str(3:end));
        IM(i).time = cell2mat(str(2));
        IM(i).values = cell2mat(str(3:end));
        
    elseif startsWith(tline,'RES') == true
        r = r + 1;
        RES(r,:) = textscan(tline,'%s %s');
        
    elseif startsWith(tline,'LANDMARK') == true
        l = l + 1;
        str = textscan(tline,'%s %f %f %f %s %d');
        LANDMARK(l).nr = cell2mat(str(6));
        LANDMARK(l).values = cell2mat(str(2:4));
        
    elseif startsWith(tline,'OBJECT') == true
        o = o + 1;
        str = textscan(tline,'%s %f %f %f %f %f %f %s %d');
        nr = cell2mat(str(9));
        OBJECT(o).nr     = cell2mat(str(9));
        OBJECT(o).values = cell2mat(str(2:7));
    end
    tline = fgetl(fid);
end
fclose(fid);

RES = cell2struct(RES(:,2),cellfun(@(x) x{1},RES(:,1),'UniformOutput',false),1);

%% Read the pixel data from the sxi file.
fid = fopen(strrep(sw_filename,'sw','sxi'));
A   = fread(fid);
fclose(fid);
PXDATA = int16(reshape(A,...
    [str2double(RES.RES_BUF_WIDTH) ...
     str2double(RES.RES_BUF_HEIGHT) ...
     str2double(RES.RES_BUF_FRAMES)]));
end