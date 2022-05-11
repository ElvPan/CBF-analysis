clear all
%%  reading data
list={'multiple tifs','single tif file','czi','lsm','lif','nd2','avi','mp4','mov','m4v'};
[indx,tf] = listdlg('PromptString','Select a file type to load','SelectionMode','single','ListString',list);
fileid=list{indx};
%reading file
switch fileid
    case {'czi','lsm','nd2','lif'}
        [FileName,PathName,FilterIndex] = uigetfile('*.czi;*.nd2;*lif;*.lsm','Select czi, lsm, nd2 or lif file');
        [image_data,specs]=readOMEFile_new(PathName,FileName);
        timeframe=specs.timeframe;
    case {'mp4','mov','m4v','avi'}
        [FileName,PathName,FilterIndex] = uigetfile('*.mp4;*.mov;*.m4v;*.avi','Select mp4,mov or m4v file');
        vidObj = VideoReader([PathName FileName]);
        vidHeight = vidObj.Height;
        vidWidth = vidObj.Width;
        s = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'));
        k = 1;
        while hasFrame(vidObj)
         s.cdata=readFrame(vidObj);
         image_data(:,:,k) = rgb2gray(s.cdata);
         k = k+1;
        end
        timeframe=[];
    case 'multiple tifs'
        PathName=uigetdir([],'Please select folder containing multiple tifs');
        S = dir(fullfile(PathName,'*'));
        N = natsortfiles({S.name});
        for k = 1:numel(N)
        image_data(:,:,k)=imread(fullfile(PathName,N{k}));
        end
        timeframe=[];
    case 'single tif file'
        [FileName,PathName,FilterIndex] = uigetfile('*.tif;*.tiff','Select multi-tif file');
        image_data=rd_img16([PathName FileName]);
        timeframe=[];
    otherwise
        warning('The file choice must be one of the available choices')
end

%% prompt which sampling approach to use to select areas of data to analyze
list={'multiple points','line','free hand','polygon'};
[indx,tf] = listdlg('PromptString','Select a sampling approach','SelectionMode','single','ListString',list);
sampleid=list{indx};
 set(gcf,'Units','normalized','Position',[.2 .2 .6 .6])
switch sampleid
    case 'multiple points'
        answer = inputdlg('How many points do you want to pick?',...
             'Sampling pts', [1 50]);
         npts=str2double(answer);
        imagesc(mean(image_data,3));colormap(gray);axis equal;axis off
        ltp=npts;
        k=1;
        clear pos
        while ltp>0
        pos1=ginput(1);
        pos(k,:)=round(pos1);
        k=k+1;
        ltp=ltp-1;
        hold on
        plot(pos1(1,1),pos1(1,2),'*g')
        end
        row=pos(:,2);
        col=pos(:,1);
        col=round(col);
        row=round(row);
        hold on
        for k=1:length(row)
           text(col(k)+10,row(k),num2str(k),'Color',[0 1 0])
           hold on
        end
    case 'line'
        imagesc(mean(image_data,3));colormap(gray);axis equal;axis off
        h=imline;
        pos = getPosition(h);
        pos=round(pos);
        ind1=find(pos(:,1)==min(pos(:,1)));
        ind2=find(pos(:,1)==max(pos(:,1)));
        x1=pos(ind1,1);
        y1=pos(ind1,2);
        x2=pos(ind2,1);
        y2=pos(ind2,2);
        a=(y2-y1)/(x2-x1);
        b=(-y2*x1+y1*x2)/(x2-x1);
        col=x1:x2;
        row=a*col+b;
        col=round(col);
        row=round(row);
        hold on
        plot(col,row,'-b')
        hold on
        for k=1:round(length(row)/10):length(row)
           plot(col(k),row(k),'*g')
           hold on
           text(col(k)+10,row(k),num2str(k),'Color',[0 1 0])
           hold on
        end
    case 'free hand'
        imagesc(mean(image_data,3));colormap(gray);axis equal;axis off
        h=imfreehand;
        pos = getPosition(h);
        row=pos(:,2);
        col=pos(:,1);
        col=round(col);
        row=round(row);
        hold on
        plot(col,row,'-b')
        hold on
        for k=1:round(length(row)/10):length(row)
           plot(col(k),row(k),'*g')
           hold on
           text(col(k)+10,row(k),num2str(k),'Color',[0 1 0])
           hold on
        end
    case 'polygon'
        imagesc(mean(image_data,3));colormap(gray);axis equal;axis off
        h=impoly;
        pos = getPosition(h);
        BW = createMask(h);
        %pos=round(pos);
        [row,col]=find(BW);
        col=round(col);
        row=round(row);
        hold on
        pos(1+size(pos,1),:)=pos(1,:);
        plot(pos(:,1),pos(:,2),'-g')
       
end
set(gcf,'Color',[1 1 1])
saveas(1,[PathName 'select_' FileName(1:end-4) '_' sampleid '.png'])
close 
%%
movie = immfilter(image_data,'F');
movie=single(movie);
%% prompt for time frame, cutoff for S/N and peak prominence 
prompt = {'Enter frame time:','Enter S/N cutoff for peak detection:','Enter peak prominence cutoff:'};
title = 'Input';
dims = [1 35];
definput = {num2str(timeframe),'5','0.5'};
answer = inputdlg(prompt,title,dims,definput);
timeframe=str2double(answer{1});
cutStoN=str2double(answer{2}); %cutoff signal to noise to detect peaks in spectrum
promCut=str2double(answer{3}); %peak prominence cutoff
%% do analysis
T = timeframe;             % Sampling period    
Fs = 1/T;            % Sampling frequency                    
L =  size(movie,3);              % Length of signal
t = (0:L-1)*T;        % Time vector
m=0;
clear modemap ampmode Peaks Locations peaksall locall specAmps specFreq maskSpectrum
%  set(gcbf,'pointer','watch');
% h = waitbar(0,'Calculating spectra per pixel and finding modes...');

for j=1:length(row)
        warning('off')
X=squeeze(movie(row(j),col(j),:));
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
 

clear movie;clear image_data
save([PathName 'workspace_' FileName(1:end-4) '.mat']);
