function [MI] = simple_motion_index(video, Timestamp, varargin)
% determine a motion index based on Jelitai et al., Nature communications,
% 2016.
% Note. HG 1/3/2017: Adapted from Fred's code.
%
% the video is a tif file (grayscale). The time stamps of the video are used to
% determine the number of frame and are appended to the matrix MI. Timestamp
% should be a one column matrix. 
% example: [MI] = Motion_index('example.tif', time).
%
% Harsha Gurnani. March 2017

MI = zeros(size(Timestamp,1),2);

tic
%% If given vector
if isnumeric(video)
        
        s=size(video);
        if size(s,2) == 3
            FrameDiff = video(:,:,2:end)-video(:,:,1:end-1);
            MI(2:end,1) = squeeze(sum(sum((FrameDiff.*FrameDiff),2),1));%^2;
        elseif size(s,2) == 2
           FrameDiff = video(:,2:end)-video(:,1:end-1);
            MI(2:end,1) = squeeze(sum((FrameDiff.*FrameDiff),1));%^2;
            
        elseif size(s,2) == 1
           FrameDiff = video(2:end)-video(1:end-1);
           MI(2:end,1) = squeeze(sum((FrameDiff.*FrameDiff),1));%^2;
        end

%If given file name
elseif ischar(video)
%% Check for substacks
    if nargin == 3
        imj_macro=0;
        n_stacks = varargin{1};

        %%% Parse path and filename
        last_bs = sort(cat(2,find(video== '/'), find(video=='\')));
        if ~isempty(last_bs)
            fpath = video(1:last_bs(end));
            flnm = video(last_bs(end)+1:end);
            if strcmp(flnm(end-3:end), '.tif')
                fltype = 'tif';
                flnm=flnm(1:end-4);
            end
            if strcmp(flnm(end-3:end), '.gif')
                fltype = 'gif';
                flnm=flnm(1:end-4);
            end
        else
            fpath = '';
            flnm = video(1:end);
        end

        if n_stacks == 0
            %%% Used IMJ macro
            imj_macro = 1;
            fls = dir([fpath,flnm, sprintf('-stack*.tiff') ]);
            n_stacks = size(fls,1);
    %         new_flnm = cell(1, n_stacks);
            for jj=1:n_stacks
            new_flnm{jj} = [fpath, flnm, sprintf('-stack%d.tiff', jj)];
    %         new_flnm{jj} = [fpath, 'stk_', sprintf('%04d',jj),'_' flnm];
            end 
        else
            new_flnm = cell(1, n_stacks);

            %%% Check filenames with '_'
            fls = dir([fpath, sprintf('*_*%s*.%s',flnm,fltype)]);
            if isempty(fls)
                fls = dir([fpath, sprintf('*%s*_*.%s',flnm,fltype)]);
            end

            for jj=1:n_stacks
                new_flnm{jj} = [fpath, fls(jj).name];
    %         new_flnm{jj} = [fpath, 'stk_', sprintf('%04d',jj),'_' flnm];
            end

        end
        max_ts = 600; %ceil(size(Timestamp,1)/n_stacks);

    else 
        n_stacks = 1;
    end

%% Compute motion index
    if n_stacks > 1

        partial_mi = nan(max_ts, n_stacks);
        parfor jj = 1:n_stacks-1
           partial_mi(:,jj) = get_part_MI( new_flnm{jj}, 1, max_ts); 
           first_frame{jj} = imread(new_flnm{jj},1);
           last_frame{jj} = imread(new_flnm{jj},max_ts);
        end

        if imj_macro
            last_id = mod(size(Timestamp,1), max_ts);
            if last_id == 0
                last_id = max_ts;
            end
            partial_mi(1:last_id,n_stacks) = get_part_MI( new_flnm{n_stacks}, 1, last_id);
            first_frame{n_stacks} = imread(new_flnm{n_stacks},1);
            last_frame{n_stacks} = imread(new_flnm{n_stacks},last_id);
        else
            partial_mi(:,n_stacks) = get_part_MI( new_flnm{n_stacks}, 1, max_ts); 
            first_frame{n_stacks} = imread(new_flnm{n_stacks},1);
            last_frame{n_stacks} = imread(new_flnm{n_stacks},max_ts);
        end

        partial_mi = reshape(partial_mi, max_ts*n_stacks,1);

        %%% Compute MI at junctions of stacks
        for jj=2:n_stacks
            FrameDiff = squeeze(first_frame{jj}(:,:,1)-last_frame{jj-1}(:,:,1));
           partial_mi(max_ts*(jj-1)+1,1) = (sum(sum(FrameDiff.*FrameDiff)));
        end
        MI(:,1) = partial_mi(1:size(Timestamp,1));

    else
        fltype = video(end-2:end);
        switch fltype
            case 'tif'
                %%% Normal case with single tif file
                parfor i = 2:size(Timestamp,1)

                    X = imread(video,i-1);
                    X1 = imread(video, i);

                    FrameDiff = squeeze(X1(:,:,1)-X(:,:,1));

                    MI(i,1) = (sum(sum(FrameDiff.*FrameDiff)));%^2;

                end

            case 'gif'
                %%% Too slow for loading in memory. DO NOT USE
                nframes = 100;
                start=1;
                ncycles = ceil(size(Timestamp,1)/nframes);
                for cyc=1:ncycles
                    last = min(start+nframes, size(Timestamp,1));
                    curr_stack = start:last;
                    vid_part = imread(video, curr_stack);
                    if cyc==1
                        temp = nan(nframes,1);
                        parfor i=2:nframes

                            FrameDiff = squeeze(vid_part(:,:,1,i-1)-vid_part(:,:,1,i));

                            temp(i,1) = (sum(sum(FrameDiff.*FrameDiff)));%^2;

                        end
                    else  
                        temp = nan(nframes,1);
                        parfor i=2:nframes+1

                            FrameDiff = squeeze(vid_part(:,:,1,i-1)-vid_part(:,:,1,i));

                            temp(i-1,1) = (sum(sum(FrameDiff.*FrameDiff)));%^2;

                        end

                    end
                    MI(curr_stack(1:end-1),1) = temp;
                    start=last;
                end

        end
    end

    %% Normalize

     m = 0;%min(MI(1:end-1000,1));
     M = max(MI(:,1));

     MI(:,1) = (MI(:,1)-m) ./ (M-m); 
     
    toc

end

MI(:,2) = Timestamp(:,2);
MI = MI(2:end,:);

end

function part_mi = get_part_MI( filenm, start, last)

   part_mi = nan(last-start+1,1);
   for i = start+1:last

        X = imread(filenm,i-1);
        X1 = imread(filenm, i);

        FrameDiff = squeeze(X1(:,:,1)-X(:,:,1));   
        part_mi(i) = (sum(sum(FrameDiff.*FrameDiff)));%^2;

    end

end