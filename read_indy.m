

function session = read_indy(session)

    StimData = zeros(session.trials_number, 8);
    tmpStimData = NaN(session.trials_number, 9+1);
    tmpDescripMatr = session.descr_mat;
    
    [fbtime, ErrorMessage] = fopen(fullfile(session.path, [session.file, '.', session.extension,'s']), 'r');
    if(fbtime<2)
	    error(['Trouble reading file ...']);
    end
    
    CInfo = fscanf(fbtime, '%s', 1);
    MInfo = fscanf(fbtime, '%s', 1);
    PInfo = fscanf(fbtime, '%s', 1);
    SInfo = fscanf(fbtime, '%s', 1);
    OInfo = fscanf(fbtime, '%s', 1);

    % Skip header line
    fgetl(fbtime);
    current_trial = 0;
    while 1
        res = fgetl(fbtime);
        
        if res==-1
            break;
        else
            res = sscanf(res, '%d\t');
        end

        if res(1) && (length(res)==res(2)+2)
            current_trial = current_trial+1;
            tmpStimData(current_trial, 1:length(res)) = res;
        end
    end
    fclose(fbtime);
    
    session.stim_data = tmpStimData;

end


