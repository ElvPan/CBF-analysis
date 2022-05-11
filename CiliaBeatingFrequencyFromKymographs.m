clear all
%%  reading data
[FileName,PathName,FilterIndex] = uigetfile('*.tif;*.tiff','Select tif files to batch process','MultiSelect','on');

if length(FileName)==1
%% read file
data=imread([PathName FileName]);
data=rgb2gray(data);
%% prompt for time frame, cutoff for S/N and peak prominence 
prompt = {'Enter frame rate (fps)','Enter S/N cutoff for peak detection:','Enter peak prominence cutoff:','Along which dimension is axis of time in images? (1=row, 2=columns)'};
title = 'Input';
dims = [1 35];
definput = {FileName(end-6:end-4),'5','0.5','2'};
answer = inputdlg(prompt,title,dims,definput);
fps=str2double(answer{1});
cutStoN=str2double(answer{2}); %cutoff signal to noise to detect peaks in spectrum
promCut=str2double(answer{3}); %peak prominence cutoff
timeaxis=answer{4};
%% see if image needs transpose to have axis of time along clumns (dimension 2)...
if strcmp(timeaxis,1)
    data=data';
end
%% do analysis
T = 1/fps;             % Sampling period    
Fs = fps;            % Sampling frequency                    
L =  size(data,2);              % Length of signal
t = (0:L-1)*T;        % Time vector
m=0;
clear modemap ampmode Peaks Locations peaksall locall specAmps specFreq maskSpectrum
%  set(gcbf,'pointer','watch');
% h = waitbar(0,'Calculating spectra per pixel and finding modes...');

for j=1:size(data,1)
        warning('off')
X=squeeze(data(j,:));
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
% 
spectouse=P1(3:end-1);
spectouse=spectouse-mean(spectouse(end-10:end-2));
freqtouse=f(3:end-1);
minHeight=median(spectouse)+2*std(spectouse);
[peak,loc]=findpeaks(spectouse,freqtouse,'MinPeakHeight',minHeight,'MinPeakDistance',3);
if ~isempty(peak)
indfilt=find(spectouse<min(peak));
minHeight=median(spectouse(indfilt))+2*std(spectouse(indfilt));
[peak,loc,w,p]=findpeaks(spectouse,freqtouse,'MinPeakHeight',minHeight,'MinPeakDistance',3);
ind1=find((peak/std(spectouse(indfilt)))>cutStoN);
ind2=find(p>promCut);
indp=intersect(ind1,ind2);
% minHeight=median(P1(2:end))+2.5*std(P1(2:end));
% [peak,loc]=findpeaks(P1(2:end),f(2:end),'MinPeakHeight',minHeight,'MinPeakDistance',1.5);
peak=peak(indp);
loc=loc(indp);
Peaks{j}=peak;
Locations{j}=loc;
peaksall(1+m:m+length(peak))=peak;
locall(1+m:m+length(peak))=loc;
specAmps(j,:)=spectouse;
specFreq(j,:)=freqtouse;

        [peaks,ind]=sort(peak);
        for k=1:length(loc)
        curfreq=loc(ind(length(loc)+1-k));
        modemap(j,k)=curfreq;
        curamp=peak(ind(length(loc)+1-k));
        ampmode(j,k)=curamp;
        end
       


        m=m+1;
%            if ishandle(h)
%                waitbar((m)/(length(row)),h)
%            else 
%                break
%            end
else
   Peaks{j}=NaN;
   Locations{j}=NaN;
   peaksall(1+m:m+length(peak))=NaN;
   locall(1+m:m+length(peak))=NaN;
   modemap(j,:)=NaN;
   ampmode(j,:)=NaN;
   specAmps(j,:)=NaN;
   specFreq(j,:)=NaN;
end

% if ishandle(h)
% close(h)
% end
% set(gcbf,'pointer','arrow');
% 
% % output average spectrum into excel and he modes extacted
% close all
end

% output excel files
clear T1
T1=table(specFreq);
writetable(T1,[PathName 'Modes_' FileName(1:end-4)  '.xlsx'],'Sheet','SpecFreq');
clear T1
T1=table(specAmps);
writetable(T1,[PathName 'Modes_' FileName(1:end-4)  '.xlsx'],'Sheet','SpecAmps');
clear T1
T1=table(modemap);
writetable(T1,[PathName 'Modes_' FileName(1:end-4)  '.xlsx'],'Sheet','fModes');
clear T
T1=table(ampmode);
writetable(T1,[PathName 'Modes_' FileName(1:end-4)  '.xlsx'],'Sheet','ampModes');
 

clear movie;clear data
save([PathName 'workspace_' FileName(1:end-4) '.mat']);
elseif length(FileName)>1
    
%% prompt for time frame, cutoff for S/N and peak prominence 
prompt = {'Enter frame rate (fps)','Enter S/N cutoff for peak detection:','Enter peak prominence cutoff:','Along which dimension is axis of time in images? (1=row, 2=columns)'};
title = 'Input';
dims = [1 35];
definput = {FileName{1}(end-6:end-4),'5','0.5','2'};
answer = inputdlg(prompt,title,dims,definput);
fps=str2double(answer{1});
cutStoN=str2double(answer{2}); %cutoff signal to noise to detect peaks in spectrum
promCut=str2double(answer{3}); %peak prominence cutoff
timeaxis=answer{4};
%% loop over files
for file=1:length(FileName)
%% read file
data=imread([PathName FileName{file}]);
data=rgb2gray(data);
%% see if image needs transpose to have axis of time along clumns (dimension 2)...
if strcmp(timeaxis,1)
    data=data';
end
%% do analysis
T = 1/fps;             % Sampling period    
Fs = fps;            % Sampling frequency                    
L =  size(data,2);              % Length of signal
t = (0:L-1)*T;        % Time vector
m=0;
clear modemap ampmode Peaks Locations peaksall locall specAmps specFreq maskSpectrum
%  set(gcbf,'pointer','watch');
% h = waitbar(0,'Calculating spectra per pixel and finding modes...');

X=squeeze(data(1,:));
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
spectouse=P1(3:end-1);

for j=1:size(data,1)
        warning('off')
X=squeeze(data(j,:));
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
% 
spectouse=P1(3:end-1);
spectouse=spectouse-mean(spectouse(end-10:end-2));
freqtouse=f(3:end-1);
minHeight=median(spectouse)+2*std(spectouse);
[peak,loc]=findpeaks(spectouse,freqtouse,'MinPeakHeight',minHeight,'MinPeakDistance',3);
if ~isempty(peak)
indfilt=find(spectouse<min(peak));
minHeight=median(spectouse(indfilt))+2*std(spectouse(indfilt));
[peak,loc,w,p]=findpeaks(spectouse,freqtouse,'MinPeakHeight',minHeight,'MinPeakDistance',3);
ind1=find((peak/std(spectouse(indfilt)))>cutStoN);
ind2=find(p>promCut);
indp=intersect(ind1,ind2);
% minHeight=median(P1(2:end))+2.5*std(P1(2:end));
% [peak,loc]=findpeaks(P1(2:end),f(2:end),'MinPeakHeight',minHeight,'MinPeakDistance',1.5);
peak=peak(indp);
loc=loc(indp);
Peaks{j}=peak;
Locations{j}=loc;
peaksall(1+m:m+length(peak))=peak;
locall(1+m:m+length(peak))=loc;
specAmps(j,:)=spectouse;
specFreq(j,:)=freqtouse;

        [peaks,ind]=sort(peak);
        for k=1:length(loc)
        curfreq=loc(ind(length(loc)+1-k));
        modemap(j,k)=curfreq;
        curamp=peak(ind(length(loc)+1-k));
        ampmode(j,k)=curamp;
        end
       


        m=m+1;
%            if ishandle(h)
%                waitbar((m)/(length(row)),h)
%            else 
%                break
%            end
else
   Peaks{j}=NaN;
   Locations{j}=NaN;
   peaksall(1+m:m+length(peak))=NaN;
   locall(1+m:m+length(peak))=NaN;
   modemap(j,:)=NaN;
   ampmode(j,:)=NaN;
   specAmps(j,:)=NaN*ones(1,length(spectouse));
   specFreq(j,:)=NaN*ones(1,length(spectouse));
end

% if ishandle(h)
% close(h)
% end
% set(gcbf,'pointer','arrow');
% 
% % output average spectrum into excel and he modes extacted
% close all
end

% output excel files
clear T1
T1=table(specFreq);
writetable(T1,[PathName 'Modes_' FileName{file}(1:end-4)  '.xlsx'],'Sheet','SpecFreq');
clear T1
T1=table(specAmps);
writetable(T1,[PathName 'Modes_' FileName{file}(1:end-4)  '.xlsx'],'Sheet','SpecAmps');
clear T1
T1=table(modemap);
writetable(T1,[PathName 'Modes_batch.xlsx'],'Sheet',['fModes' FileName{file}(12:end-4)]);
clear T
T1=table(ampmode);
writetable(T1,[PathName 'Modes_batch.xlsx'],'Sheet',['ampModes' FileName{file}(12:end-4)]);
 

clear movie;clear data
save([PathName 'workspace_' FileName{file}(1:end-4) '.mat']);  
end
end