
function session = read_data(session)
    

    [fbdata, ErrorMessage] = fopen(fullfile(session.path, [session.file, '.', session.extension,'d']), 'rb');
    if(fbdata<2)
	    error(ErrorMeassage);
    end
    
    fseek(fbdata, 0, 'eof');
    FileSize = ftell(fbdata);
    fseek(fbdata, 0, 'bof');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
    %%%%%%%%	Read Header		%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
    
    session.date = fread(fbdata, 3, 'short');
    session.time = fread(fbdata, 2, 'short');
    session.comment = fread(fbdata, 80, 'char');
    session.numEcode = fread(fbdata, 1, 'short');
    MaxArrayCode = session.numEcode;
    
    session.trials_number = fread(fbdata, 1, 'short');
    %Achtung: Die folgende Zeile war nur eingefügt, weil bei einer Aufnahme der Pen025
    %ein Mismatch zwischen INDY file und Nabeda Daten nicht anders zu beheben war. 
    %NumTrial = 209;
    
    HeaderDummy = fread(fbdata, 1, 'short');
    session.spike_channels_number = fread(fbdata, 1, 'short');   
    session.spike_channels_specs = fread(fbdata, 5*session.spike_channels_number, 'short');
%     set(ck_Sp1, 'value', 1); set(ck_Sp2, 'value', 1); set(ck_Sp3, 'value', 1); set(ck_Sp4, 'value', 1);
% 	set(ck_Sp5, 'value', 1); set(ck_Sp6, 'value', 1); set(ck_Sp7, 'value', 1); set(ck_Sp8, 'value', 1);
    
    session.analog_channels_number = fread(fbdata, 1, 'short');
    
    if(session.analog_channels_number>4)
	    error('This program is dedicated only for data with up to four analog channels, sorry. Modify READDATA.M ...');
    end
    
    session.analog_channels_samplerate = NaN(1, session.analog_channels_number);
    session.analog_channels_session.analog_channels_gain = NaN(1, session.analog_channels_number);
    for i = 1:session.analog_channels_number
	    session.analog_channels_samplerate(i) = fread(fbdata, 1, 'short');
	    session.analog_channels_gain(i) = fread(fbdata, 1, 'float');
    end
    
    session.header_size = ftell(fbdata);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
    %%%%%%%%	Read Trial by Trial		%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
    
    NumGlobalParam = zeros(session.trials_number,1);
    NumEvent = zeros(session.trials_number,1);
    NumSpike = zeros(session.trials_number,session.spike_channels_number);
    NumAnaData = zeros(session.trials_number,session.analog_channels_number);
    AnaRange2 = zeros(session.trials_number,session.analog_channels_number);
    FPosSpDat1 = zeros(session.trials_number,1);
    FPosSpDat2 = zeros(session.trials_number,1);
    
    MaxNumSpike = zeros(1,session.spike_channels_number);
    MaxAnaData = zeros(1,session.analog_channels_number);
    
    for tmpNumTrial = 1:session.trials_number
	    NumGlobalParam(tmpNumTrial,1) = fread(fbdata, 1, 'short');
	    Dummy = fread(fbdata, NumGlobalParam(tmpNumTrial,1), 'short');
	    NumEvent(tmpNumTrial,1) = fread(fbdata, 1, 'short');
	    for tmpNumEvent = 1:NumEvent(tmpNumTrial,1)
		    Dummy = fread(fbdata, 1, 'ulong');
		    NumEvParam = fread(fbdata, 1, 'short');
		    EvParam = fread(fbdata, NumEvParam, 'short');
        end
	    for tmpNumSpikeChan = 1:session.spike_channels_number
		    NumSpike(tmpNumTrial, tmpNumSpikeChan) = fread(fbdata, 1, 'ulong');
		    Dummy = fread(fbdata, NumSpike(tmpNumTrial, tmpNumSpikeChan), 'short');
        end
	    FPosSpDat1(tmpNumTrial) = ftell(fbdata);
	    for tmpNumAnaChan = 1:session.analog_channels_number
		    NumAnaData(tmpNumTrial, tmpNumAnaChan) = fread(fbdata, 1, 'ulong');
		    Dummy = fread(fbdata, NumAnaData(tmpNumTrial, tmpNumAnaChan), 'short');
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ReadSize = ftell(fbdata);
    
    if(ReadSize ~= FileSize)
	    error('File size mismatch ...');
    end
    
    
    for index = 1:session.spike_channels_number
	    MaxNumSpike(1,index) = max(NumSpike(:,index));
    end
    for index = 1:session.analog_channels_number
	    MaxAnaData(index) = max(NumAnaData(:,index));
    end
    
    GlobalParam = zeros(session.trials_number, max(max(NumGlobalParam,1)));
    EventTime = zeros(session.trials_number, max(max(NumEvent,1)));
    if(session.spike_channels_number==1)
	    Spike1Time = zeros(session.trials_number, MaxNumSpike(1,1));
    elseif(session.spike_channels_number==2)
	    Spike1Time = zeros(session.trials_number, MaxNumSpike(1,1));
	    Spike2Time = zeros(session.trials_number, MaxNumSpike(1,2));
    elseif(session.spike_channels_number==3)
	    Spike1Time = zeros(session.trials_number, MaxNumSpike(1,1));
	    Spike2Time = zeros(session.trials_number, MaxNumSpike(1,2));
	    Spike3Time = zeros(session.trials_number, MaxNumSpike(1,3));
    elseif(session.spike_channels_number==4)
	    Spike1Time = zeros(session.trials_number, MaxNumSpike(1,1));
	    Spike2Time = zeros(session.trials_number, MaxNumSpike(1,2));
	    Spike3Time = zeros(session.trials_number, MaxNumSpike(1,3));
	    Spike4Time = zeros(session.trials_number, MaxNumSpike(1,4));
    elseif(session.spike_channels_number==5)
	    Spike1Time = zeros(session.trials_number, MaxNumSpike(1,1));
	    Spike2Time = zeros(session.trials_number, MaxNumSpike(1,2));
	    Spike3Time = zeros(session.trials_number, MaxNumSpike(1,3));
	    Spike4Time = zeros(session.trials_number, MaxNumSpike(1,4));
	    Spike5Time = zeros(session.trials_number, MaxNumSpike(1,5));
    elseif(session.spike_channels_number==6)
	    Spike1Time = zeros(session.trials_number, MaxNumSpike(1,1));
	    Spike2Time = zeros(session.trials_number, MaxNumSpike(1,2));
	    Spike3Time = zeros(session.trials_number, MaxNumSpike(1,3));
	    Spike4Time = zeros(session.trials_number, MaxNumSpike(1,4));
	    Spike5Time = zeros(session.trials_number, MaxNumSpike(1,5));
	    Spike6Time = zeros(session.trials_number, MaxNumSpike(1,6));
    elseif(session.spike_channels_number==7)
	    Spike1Time = zeros(session.trials_number, MaxNumSpike(1,1));
	    Spike2Time = zeros(session.trials_number, MaxNumSpike(1,2));
	    Spike3Time = zeros(session.trials_number, MaxNumSpike(1,3));
	    Spike4Time = zeros(session.trials_number, MaxNumSpike(1,4));
	    Spike5Time = zeros(session.trials_number, MaxNumSpike(1,5));
	    Spike6Time = zeros(session.trials_number, MaxNumSpike(1,6));
	    Spike7Time = zeros(session.trials_number, MaxNumSpike(1,7));
    elseif(session.spike_channels_number>7)
	    error('This program is dedicated only for data with up to seven spike channels, sorry ...');
    end
    
    BoolAna(1) = 1;
    BoolAna(2) = 1;
    BoolAna(3) = 1;
    BoolAna(4) = 1;
    Check = sum(BoolAna);
    if(Check>0)
	    if((BoolAna(1)==1))
		    Ana1Data = zeros(session.trials_number, MaxAnaData(1,1));
        end
	    if((BoolAna(2)==1))
		    Ana2Data = zeros(session.trials_number, MaxAnaData(1,2));
        end
	    if((BoolAna(3)==1))
		    Ana3Data = zeros(session.trials_number, MaxAnaData(1,3));
        end
	    if((BoolAna(4)==1))
		    Ana4Data = zeros(session.trials_number, MaxAnaData(1,4));
        end
    end
    
    fseek(fbdata, session.header_size, 'bof');
    FilePos = ftell(fbdata);
    if(FilePos ~= session.header_size)
	    error('File read mismatch ...');
    end
    for tmpNumTrial = 1:session.trials_number
	    DummyRange = fread(fbdata, 1, 'short');
	    GlobalParam(tmpNumTrial,1:NumGlobalParam(tmpNumTrial,1)) = fread(fbdata, NumGlobalParam(tmpNumTrial,1), 'short')';
	    DummyRange = fread(fbdata, 1, 'short');
	    for tmpNumEvent = 1:NumEvent(tmpNumTrial,1)
		    EventTime(tmpNumTrial,tmpNumEvent) = fread(fbdata, 1, 'ulong');
		    Dummy1 = fread(fbdata, 1, 'short');
		    Dummy2 = fread(fbdata, Dummy1, 'short');
        end
	    for tmpNumSpikeChan = 1:session.spike_channels_number
		    Range = fread(fbdata, 1, 'ulong');
		    if(Range ~= NumSpike(tmpNumTrial, tmpNumSpikeChan))
			    error('NumSpike mismatch ...');
            end
		    if(tmpNumSpikeChan==1)
			    Spike1Time(tmpNumTrial,1:Range) = fread(fbdata, Range, 'short')'*(1000.0/session.spike_channels_specs(1));
		    elseif(tmpNumSpikeChan==2)
			    Spike2Time(tmpNumTrial,1:Range) = fread(fbdata, Range, 'short')'*(1000.0/session.spike_channels_specs(6));
		    elseif(tmpNumSpikeChan==3)
			    Spike3Time(tmpNumTrial,1:Range) = fread(fbdata, Range, 'short')'*(1000.0/session.spike_channels_specs(11));
		    elseif(tmpNumSpikeChan==4)
			    Spike4Time(tmpNumTrial,1:Range) = fread(fbdata, Range, 'short')'*(1000.0/session.spike_channels_specs(16));
		    elseif(tmpNumSpikeChan==5)
			    Spike5Time(tmpNumTrial,1:Range) = fread(fbdata, Range, 'short')'*(1000.0/session.spike_channels_specs(21));
		    elseif(tmpNumSpikeChan==6)
			    Spike6Time(tmpNumTrial,1:Range) = fread(fbdata, Range, 'short')'*(1000.0/session.spike_channels_specs(26));
		    elseif(tmpNumSpikeChan==7)
			    Spike7Time(tmpNumTrial,1:Range) = fread(fbdata, Range, 'short')'*(1000.0/session.spike_channels_specs(31));
            end
        end
    
	    FPosSpDat2(tmpNumTrial) = ftell(fbdata);
	    for tmpNumAnaChan = 1:session.analog_channels_number
		    AnaRange2(tmpNumTrial,tmpNumAnaChan) = fread(fbdata, 1, 'ulong');
		    Range = AnaRange2(tmpNumTrial,tmpNumAnaChan);
		    if(Range ~= NumAnaData(tmpNumTrial, tmpNumAnaChan))
			    error('NumAnaData mismatch ...');
            end
		    if(tmpNumAnaChan==1) 
			    if(BoolAna(1)==1)
				    Ana1Data(tmpNumTrial,1:Range) = fread(fbdata, Range, 'short')'/session.analog_channels_gain(1);
			    else
				    Dummy = fread(fbdata, Range, 'short')';
                end
    
		    elseif(tmpNumAnaChan==2) 
			    if(BoolAna(2)==1)
				    Ana2Data(tmpNumTrial,1:Range) = fread(fbdata, Range, 'short')'/session.analog_channels_gain(2);
			    else
				    Dummy = fread(fbdata, Range, 'short')';
                end
    
		    elseif(tmpNumAnaChan==3) 
			    if(BoolAna(3)==1)
				    Ana3Data(tmpNumTrial,1:Range) = fread(fbdata, Range, 'short')'/session.analog_channels_gain(3);
			    else
				    Dummy = fread(fbdata, Range, 'short')';
                end
    
		    elseif(tmpNumAnaChan==4) 
			    if(BoolAna(4)==1)
				    Ana4Data(tmpNumTrial,1:Range) = fread(fbdata, Range, 'short')'/session.analog_channels_gain(4);
			    else
				    Dummy = fread(fbdata, Range, 'short')';
                end
            end
        end
    end
    
    max_range = max(NumAnaData,[],'all');
    session.analog_data = NaN(session.trials_number, 2*max_range, 4);
    for trial=1:session.trials_number
        session.analog_data(trial,1:2*NumAnaData(trial,1),1) = reshape(repmat(Ana1Data(trial,1:NumAnaData(trial,1)),2,1),1,[]);
        session.analog_data(trial,1:2*NumAnaData(trial,2),2) = reshape(repmat(Ana2Data(trial,1:NumAnaData(trial,2)),2,1),1,[]);
        session.analog_data(trial,1:2*NumAnaData(trial,3),3) = reshape(repmat(Ana3Data(trial,1:NumAnaData(trial,3)),2,1),1,[]);
        session.analog_data(trial,1:2*NumAnaData(trial,4),4) = reshape(repmat(Ana4Data(trial,1:NumAnaData(trial,4)),2,1),1,[]);
    end
    
    try
        session.spike_data = NaN(session.trials_number, 2*max_range, 7);
        for trial=1:session.trials_number
            session.spike_data(trial,1:2*NumAnaData(trial,1),:) = 0;
    
            if NumSpike(trial,1)>0
                session.spike_data(trial,max(Spike1Time(trial,1:NumSpike(trial,1)),1),1) = 1;
            end
            if NumSpike(trial,2)>0
                session.spike_data(trial,max(Spike2Time(trial,1:NumSpike(trial,2)),1),2) = 1;
            end
            if NumSpike(trial,3)>0
                session.spike_data(trial,max(Spike3Time(trial,1:NumSpike(trial,3)),1),3) = 1;
            end
            if NumSpike(trial,4)>0
                session.spike_data(trial,max(Spike4Time(trial,1:NumSpike(trial,4)),1),4) = 1;
            end
            if NumSpike(trial,5)>0
                session.spike_data(trial,max(Spike5Time(trial,1:NumSpike(trial,5)),1),5) = 1;
            end
            if NumSpike(trial,6)>0
                session.spike_data(trial,max(Spike6Time(trial,1:NumSpike(trial,6)),1),6) = 1;
            end
            if NumSpike(trial,7)>0
                session.spike_data(trial,max(Spike7Time(trial,1:NumSpike(trial,7)),1),7) = 1;
            end
        end
    catch
        disp('oops')
    end
    
    session.events_time = EventTime;
    session.global_params = GlobalParam;
    
    session.descr_mat = NaN(session.trials_number,17);
    session.descr_mat(:,1) = (1:session.trials_number);
    session.descr_mat(:,2) = session.global_params(:,1);
    session.descr_mat(:,3) = ones(session.trials_number,1);
    if(size(session.events_time,2)==5)
	    session.descr_mat(:,6:10) = session.events_time;
    elseif(size(session.events_time,2)==4)
	    session.descr_mat(:,6:9) = session.events_time;
    end

    fclose(fbdata);

end

