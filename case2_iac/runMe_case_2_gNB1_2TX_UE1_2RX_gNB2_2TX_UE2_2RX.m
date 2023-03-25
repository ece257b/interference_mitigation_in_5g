%% Case 2 - Interference Alignment and Cancellation
%% Simulation Parameters
% Specify the signal-to-noise ratio (SNR), number of slots to simulate, and perfect channel estimation flag

fclose all;clear;clc;close all;

snrRange = -10:20; % SNR in dB
numMC = 5;
totalNoSlots = 1; % Number of slots to simulate
rng("default");
noiseless = false; % add noise at UE
plotFigures = false;
num_gNB = 2;
num_UE  = 2;


%% Carrier Configuration
carrier = nrCarrierConfig;
carrier.NSizeGrid = 52;

simTimer = tic;
BER_IAC = zeros(length(snrRange),num_UE,totalNoSlots);
BER_INT = zeros(length(snrRange),num_UE,totalNoSlots);
BER_NO_INT = zeros(length(snrRange),num_UE,totalNoSlots);
for snrIdx = 1:length(snrRange)
  SNRdB = snrRange(snrIdx);
  berIac = zeros(num_UE,totalNoSlots);
  berInt = zeros(num_UE,totalNoSlots);
  berNoInt = zeros(num_UE,totalNoSlots);
  for mcIdx = 1:numMC

    %% PDSCH and DM-RS Configuration
    pdsch = cell(1,num_gNB);
    codeRate = cell(1,num_gNB);
    encodeDLSCH = cell(1,num_gNB);
    encodeDLSCHInt = cell(1,num_gNB);
    encodeDLSCHNoInt = cell(1,num_gNB);
    decodeDLSCH = cell(1,num_gNB);
    decodeDLSCHInt = cell(1,num_gNB);
    decodeDLSCHNoInt = cell(1,num_gNB);
    harqEntity = cell(1,num_gNB);
    harqEntityInt = cell(1,num_gNB);
    harqEntityNoInt = cell(1,num_gNB);
    for gnbIdx = 1:num_gNB
      %% PDSCH and DM-RS Configuration
      pdsch{gnbIdx} = nrPDSCHConfig;
      pdsch{gnbIdx}.Modulation = "QPSK";
      % pdsch_1.Modulation = "16QAM";
%       pdsch{gnbIdx}.NumLayers = 1;
      pdsch{gnbIdx}.PRBSet = 0:carrier.NSizeGrid-1; % Full band allocation

      % Set the DM-RS parameters. To improve channel estimation, add an additional DM-RS position
      pdsch{gnbIdx}.DMRS.DMRSAdditionalPosition = 1;

      % Set the DM-RS configuration type and the DM-RS length, which determines the number of orthogonal DM-RS sequences or DM-RS ports.
      % The maximum number of layers must be less than or equal to the number of DM-RS ports.
      pdsch{gnbIdx}.DMRS.DMRSConfigurationType = 1;
      pdsch{gnbIdx}.DMRS.DMRSLength = 2;


      %% DL-SCH Configuration
      NHARQProcesses = 16; % Number of parallel HARQ processes
      rvSeq = [0 2 3 1];

      % Coding rate
      if pdsch{gnbIdx}.NumCodewords == 1
        codeRate{gnbIdx} = 490/1024;
      else
        codeRate{gnbIdx} = [490 490]./1024;
      end

      % Create DL-SCH encoder object
      encodeDLSCH{gnbIdx} = nrDLSCH;
      encodeDLSCH{gnbIdx}.MultipleHARQProcesses = true;
      encodeDLSCH{gnbIdx}.TargetCodeRate = codeRate{gnbIdx};

      encodeDLSCHInt{gnbIdx} = nrDLSCH;
      encodeDLSCHInt{gnbIdx}.MultipleHARQProcesses = true;
      encodeDLSCHInt{gnbIdx}.TargetCodeRate = codeRate{gnbIdx};

      encodeDLSCHNoInt{gnbIdx} = nrDLSCH;
      encodeDLSCHNoInt{gnbIdx}.MultipleHARQProcesses = true;
      encodeDLSCHNoInt{gnbIdx}.TargetCodeRate = codeRate{gnbIdx};

      % Create DLSCH decoder object
      decodeDLSCH{gnbIdx} = nrDLSCHDecoder;
      decodeDLSCH{gnbIdx}.MultipleHARQProcesses = true;
      decodeDLSCH{gnbIdx}.TargetCodeRate = codeRate{gnbIdx};
      decodeDLSCH{gnbIdx}.LDPCDecodingAlgorithm = "Normalized min-sum";
      decodeDLSCH{gnbIdx}.MaximumLDPCIterationCount = 6;

      decodeDLSCHInt{gnbIdx} = nrDLSCHDecoder;
      decodeDLSCHInt{gnbIdx}.MultipleHARQProcesses = true;
      decodeDLSCHInt{gnbIdx}.TargetCodeRate = codeRate{gnbIdx};
      decodeDLSCHInt{gnbIdx}.LDPCDecodingAlgorithm = "Normalized min-sum";
      decodeDLSCHInt{gnbIdx}.MaximumLDPCIterationCount = 6;

      decodeDLSCHNoInt{gnbIdx} = nrDLSCHDecoder;
      decodeDLSCHNoInt{gnbIdx}.MultipleHARQProcesses = true;
      decodeDLSCHNoInt{gnbIdx}.TargetCodeRate = codeRate{gnbIdx};
      decodeDLSCHNoInt{gnbIdx}.LDPCDecodingAlgorithm = "Normalized min-sum";
      decodeDLSCHNoInt{gnbIdx}.MaximumLDPCIterationCount = 6;


      %% HARQ Management
      harqEntity{gnbIdx} = HARQEntity(0:NHARQProcesses-1,rvSeq,pdsch{gnbIdx}.NumCodewords);
      harqEntityInt{gnbIdx} = HARQEntity(0:NHARQProcesses-1,rvSeq,pdsch{gnbIdx}.NumCodewords);
      harqEntityNoInt{gnbIdx} = HARQEntity(0:NHARQProcesses-1,rvSeq,pdsch{gnbIdx}.NumCodewords);
    end
    pdsch{1}.NumLayers = 1;
    pdsch{2}.NumLayers = 2;


    %% Channel Configuration
    % Specify the number of transmit and receive antennas
    gNB_numTxAnts = [2, 2];
    UE_numRxAnts = [2, 2];

    ofdmInfo = nrOFDMInfo(carrier);

    % channel from gNB x to UE y
    seed = 1;
    channel = cell(num_gNB,num_UE);
    for ueIdx = 1:num_UE
      for gnbIdx = 1:num_gNB
        channel{gnbIdx,ueIdx} = nrTDLChannel;
        channel{gnbIdx,ueIdx}.DelayProfile = 'TDL-A';
        channel{gnbIdx,ueIdx}.NumTransmitAntennas = gNB_numTxAnts(gnbIdx);
        channel{gnbIdx,ueIdx}.NumReceiveAntennas = UE_numRxAnts(ueIdx);
        channel{gnbIdx,ueIdx}.SampleRate = ofdmInfo.SampleRate;
        channel{gnbIdx,ueIdx}.Seed = seed;
        seed = seed + 1;
      end
    end

    %% gNB Precoding Weight Calculations
    estChannelGrid = cell(num_gNB,num_UE);
    for gnbIdx = 1:num_gNB
      for ueIdx = 1:num_UE
        estChannelGrid{gnbIdx,ueIdx} = getInitialChannelEstimate(channel{gnbIdx,ueIdx},carrier);
      end
    end

    % obtain precoding weights for gNB1 (align x_1 with x_2 to allow decoding of x_3 at UE2)
    precodingWeights{1} = getPrecodingMatrixIA(pdsch{1}.PRBSet,estChannelGrid{1,2},estChannelGrid{2,2});

    % obtain precoding weights for gNB1 (non IAC)
    precodingWeightsIdeal{1} = [1 1] / sqrt(2);

    % obtain precoding weights for gNB2 (non IAC)
    precodingWeightsIdeal{2} = eye(pdsch{2}.NumLayers);


    %% Transmission and Reception
    % Simulate the transmission and reception of slots
    offset = 0;
    decbits = cell(totalNoSlots,num_UE);
    blkerr = cell(totalNoSlots,num_UE);
    for nSlot = 0:totalNoSlots-1
      % New slot
      carrier.NSlot = nSlot;
      % Generate PDSCH indices and info
      pdschIndices = cell(1,num_gNB);
      pdschInfo = cell(1,num_gNB);
      for gnbIdx = 1:num_gNB
        [pdschIndices{gnbIdx},pdschInfo{gnbIdx}] = nrPDSCHIndices(carrier,pdsch{gnbIdx});
      end


      %% gNB - TB generation, DL-SCH Encoding, PDSCH Modulation, and MIMO / Beam Null / IAC Precoding
      trBlkSizes = cell(1,num_gNB);
      codedTrBlock = cell(1,num_gNB);
      codedTrBlockInt = cell(1,num_gNB);
      codedTrBlockNoInt = cell(1,num_gNB);
      pdschSymbolsRaw = cell(1,num_gNB);
      pdschSymbolsRawInt = cell(1,num_gNB);
      pdschSymbolsRawNoInt = cell(1,num_gNB);
      pdschSymbols = cell(1,num_gNB);
      pdschSymbolsInt = cell(1,num_gNB);
      pdschSymbolsNoInt = cell(1,num_gNB);
      dmrsSymbols = cell(1,num_gNB);
      dmrsIndices = cell(1,num_gNB);
      pdschGrid = cell(1,num_gNB);
      pdschGridInt = cell(1,num_gNB);
      pdschGridNoInt = cell(1,num_gNB);
      pdschAntIndices = cell(1,num_gNB);
      dmrsAntIndices = cell(1,num_gNB);
      for gnbIdx = 1:num_gNB
        % Calculate transport block sizes
        Xoh_PDSCH = 0;
        trBlkSizes{gnbIdx} = nrTBS(pdsch{gnbIdx}.Modulation,pdsch{gnbIdx}.NumLayers,numel(pdsch{gnbIdx}.PRBSet),pdschInfo{gnbIdx}.NREPerPRB,codeRate{gnbIdx},Xoh_PDSCH);

        % HARQ Processing
        % Get new transport blocks and flush decoder soft buffer, as required
        for cwIdx = 1:pdsch{gnbIdx}.NumCodewords
          if harqEntity{gnbIdx}.NewData(cwIdx)
            % Create and store a new transport block for transmission
            trBlk = randi([0 1],trBlkSizes{gnbIdx}(cwIdx),1);
            setTransportBlock(encodeDLSCH{gnbIdx},trBlk,cwIdx-1,harqEntity{gnbIdx}.HARQProcessID);
            setTransportBlock(encodeDLSCHInt{gnbIdx},trBlk,cwIdx-1,harqEntityInt{gnbIdx}.HARQProcessID);
            setTransportBlock(encodeDLSCHNoInt{gnbIdx},trBlk,cwIdx-1,harqEntityNoInt{gnbIdx}.HARQProcessID);

            % If the previous RV sequence ends without successful decoding, flush the soft buffer
            if harqEntity{gnbIdx}.SequenceTimeout(cwIdx)
              resetSoftBuffer(decodeDLSCH{gnbIdx},cwIdx-1,harqEntity{gnbIdx}.HARQProcessID);
            end
            if harqEntityInt{gnbIdx}.SequenceTimeout(cwIdx)
              resetSoftBuffer(decodeDLSCHInt{gnbIdx},cwIdx-1,harqEntityInt{gnbIdx}.HARQProcessID);
            end
            if harqEntityNoInt{gnbIdx}.SequenceTimeout(cwIdx)
              resetSoftBuffer(decodeDLSCHNoInt{gnbIdx},cwIdx-1,harqEntityNoInt{gnbIdx}.HARQProcessID);
            end
          end
        end

        % DL-SCH Encoding
        codedTrBlock{gnbIdx} = encodeDLSCH{gnbIdx}(pdsch{gnbIdx}.Modulation,pdsch{gnbIdx}.NumLayers,pdschInfo{gnbIdx}.G,harqEntity{gnbIdx}.RedundancyVersion,harqEntity{gnbIdx}.HARQProcessID);
        codedTrBlockInt{gnbIdx} = encodeDLSCHInt{gnbIdx}(pdsch{gnbIdx}.Modulation,pdsch{gnbIdx}.NumLayers,pdschInfo{gnbIdx}.G,harqEntityInt{gnbIdx}.RedundancyVersion,harqEntityInt{gnbIdx}.HARQProcessID);
        codedTrBlockNoInt{gnbIdx} = encodeDLSCHNoInt{gnbIdx}(pdsch{gnbIdx}.Modulation,pdsch{gnbIdx}.NumLayers,pdschInfo{gnbIdx}.G,harqEntityNoInt{gnbIdx}.RedundancyVersion,harqEntityNoInt{gnbIdx}.HARQProcessID);

        % PDSCH Modulation and Precoding
        pdschSymbolsRaw{gnbIdx} = nrPDSCH(carrier,pdsch{gnbIdx},codedTrBlock{gnbIdx}); % (xdim = populated subcarriers * 14 symbols, ydim = numLayers)
        pdschSymbolsRawInt{gnbIdx} = nrPDSCH(carrier,pdsch{gnbIdx},codedTrBlockInt{gnbIdx}); % (xdim = populated subcarriers * 14 symbols, ydim = numLayers)
        pdschSymbolsRawNoInt{gnbIdx} = nrPDSCH(carrier,pdsch{gnbIdx},codedTrBlockNoInt{gnbIdx}); % (xdim = populated subcarriers * 14 symbols, ydim = numLayers)

        switch gnbIdx
          case 1
            % gNB 1 - IA precoding
            % Precode the PDSCH symbols (beam nulling done after inserting into resource grid to allow for per-subcarrier precoding)
            pdschSymbols{gnbIdx} = pdschSymbolsRaw{gnbIdx}*[1, 1] / sqrt(2); % create 2 copies of pdschSymbols (setup for precoding)

            % PDSCH DM-RS Generation
            dmrsSymbols{gnbIdx} = nrPDSCHDMRS(carrier,pdsch{gnbIdx});
            dmrsIndices{gnbIdx} = nrPDSCHDMRSIndices(carrier,pdsch{gnbIdx});

            % gNB 1 - Mapping to Resource Grid
            pdschGrid{gnbIdx} = nrResourceGrid(carrier,gNB_numTxAnts(gnbIdx));

            % insert PDSCH symbols into PDSCH resourse grid
            [~,pdschAntIndices{gnbIdx}] = nrExtractResources(pdschIndices{gnbIdx},pdschGrid{gnbIdx});
            pdschGrid{gnbIdx}(pdschAntIndices{gnbIdx}) = pdschSymbols{gnbIdx};

            % Precode and map the DM-RS symbols to the resource grid
            for p = 1:size(dmrsSymbols{gnbIdx},2)
              [~,dmrsAntIndices{gnbIdx}] = nrExtractResources(dmrsIndices{gnbIdx}(:,p),pdschGrid{gnbIdx});
              pdschGrid{gnbIdx}(dmrsAntIndices{gnbIdx}) = pdschGrid{gnbIdx}(dmrsAntIndices{gnbIdx}) + dmrsSymbols{gnbIdx}(:,p); % dmrsSymbols will be precoded below
            end

            % gNB 1
            for symIdx = 1:size(pdschGrid{gnbIdx},2)
              [~,~,A] = size(pdschGrid{gnbIdx});
              symGrid = reshape(pdschGrid{gnbIdx}(:,symIdx,:),[],A);
              pdschGrid{gnbIdx}(:,symIdx,:) = symGrid .* precodingWeights{gnbIdx};
            end
          case 2
            % gNB 2 - No Precoding
            % Precode the PDSCH symbols (beam nulling done after inserting into resource grid to allow for per-subcarrier precoding)
            pdschSymbols{gnbIdx} = pdschSymbolsRaw{gnbIdx} / sqrt(2); % create 2 copies of pdschSymbols (setup for precoding)

            % PDSCH DM-RS Generation
            dmrsSymbols{gnbIdx} = nrPDSCHDMRS(carrier,pdsch{gnbIdx});
            dmrsIndices{gnbIdx} = nrPDSCHDMRSIndices(carrier,pdsch{gnbIdx});

            % gNB 2 - Mapping to Resource Grid
            pdschGrid{gnbIdx} = nrResourceGrid(carrier,gNB_numTxAnts(gnbIdx));

            % insert PDSCH symbols into PDSCH resourse grid
            [~,pdschAntIndices{gnbIdx}] = nrExtractResources(pdschIndices{gnbIdx},pdschGrid{gnbIdx});
            pdschGrid{gnbIdx}(pdschAntIndices{gnbIdx}) = pdschSymbols{gnbIdx};

            % Precode and map the DM-RS symbols to the resource grid
            for p = 1:size(dmrsSymbols{gnbIdx},2)
              [~,dmrsAntIndices{gnbIdx}] = nrExtractResources(dmrsIndices{gnbIdx}(:,p),pdschGrid{gnbIdx});
              pdschGrid{gnbIdx}(dmrsAntIndices{gnbIdx}) = pdschGrid{gnbIdx}(dmrsAntIndices{gnbIdx}) + dmrsSymbols{gnbIdx}(:,p);
            end
          case 3
            % gNB 3 - IAC Precoding
            % Precode the PDSCH symbols (beam nulling done after inserting into resource grid to allow for per-subcarrier precoding)
            pdschSymbols{gnbIdx} = pdschSymbolsRaw{gnbIdx}*[1, 1, 1]; % create 3 copies of pdschSymbols (setup for precoding)

            % PDSCH DM-RS Generation
            dmrsSymbols{gnbIdx} = nrPDSCHDMRS(carrier,pdsch{gnbIdx});
            dmrsIndices{gnbIdx} = nrPDSCHDMRSIndices(carrier,pdsch{gnbIdx});

            % gNB 3 - Mapping to Resource Grid
            pdschGrid{gnbIdx} = nrResourceGrid(carrier,gNB_numTxAnts(gnbIdx));

            % insert PDSCH symbols into PDSCH resourse grid
            [~,pdschAntIndices{gnbIdx}] = nrExtractResources(pdschIndices{gnbIdx},pdschGrid{gnbIdx});
            pdschGrid{gnbIdx}(pdschAntIndices{gnbIdx}) = pdschSymbols{gnbIdx};

            % Precode and map the DM-RS symbols to the resource grid
            for p = 1:size(dmrsSymbols{gnbIdx},2)
              [~,dmrsAntIndices{gnbIdx}] = nrExtractResources(dmrsIndices{gnbIdx}(:,p),pdschGrid{gnbIdx});
              pdschGrid{gnbIdx}(dmrsAntIndices{gnbIdx}) = pdschGrid{gnbIdx}(dmrsAntIndices{gnbIdx}) + dmrsSymbols{gnbIdx}(:,p); % dmrsSymbols will be precoded below
            end

            % gNB 3 IAC precoding
            for symIdx = 1:size(pdschGrid{gnbIdx},2)
              [~,~,A] = size(pdschGrid{gnbIdx});
              symGrid = reshape(pdschGrid{gnbIdx}(:,symIdx,:),[],A);
              pdschGrid{gnbIdx}(:,symIdx,:) = symGrid .* precodingWeights{gnbIdx};
            end
          otherwise
            assert(0);
        end

        %% Int
        % No precoding (antenna scaling)
        pdschSymbolsInt{gnbIdx} = pdschSymbolsRawInt{gnbIdx}*precodingWeightsIdeal{gnbIdx};

        % Mapping to Resource Grid
        pdschGridInt{gnbIdx} = nrResourceGrid(carrier,gNB_numTxAnts(gnbIdx));

        % insert PDSCH symbols into PDSCH resourse grid
        pdschGridInt{gnbIdx}(pdschAntIndices{gnbIdx}) = pdschSymbolsInt{gnbIdx};

        % Precode and map the DM-RS symbols to the resource grid
        for p = 1:size(dmrsSymbols{gnbIdx},2)
          [~,dmrsAntIndices{gnbIdx}] = nrExtractResources(dmrsIndices{gnbIdx}(:,p),pdschGridInt{gnbIdx});
          pdschGridInt{gnbIdx}(dmrsAntIndices{gnbIdx}) = pdschGridInt{gnbIdx}(dmrsAntIndices{gnbIdx}) + dmrsSymbols{gnbIdx}(:,p)*precodingWeightsIdeal{gnbIdx}(p,:);
        end

        %% No Int
        % No precoding (antenna scaling)
        pdschSymbolsNoInt{gnbIdx} = pdschSymbolsRawNoInt{gnbIdx}*precodingWeightsIdeal{gnbIdx};

        % Mapping to Resource Grid
        pdschGridNoInt{gnbIdx} = nrResourceGrid(carrier,gNB_numTxAnts(gnbIdx));

        % insert PDSCH symbols into PDSCH resourse grid
        pdschGridNoInt{gnbIdx}(pdschAntIndices{gnbIdx}) = pdschSymbolsNoInt{gnbIdx};

        % Precode and map the DM-RS symbols to the resource grid
        for p = 1:size(dmrsSymbols{gnbIdx},2)
          [~,dmrsAntIndices{gnbIdx}] = nrExtractResources(dmrsIndices{gnbIdx}(:,p),pdschGridNoInt{gnbIdx});
          pdschGridNoInt{gnbIdx}(dmrsAntIndices{gnbIdx}) = pdschGridNoInt{gnbIdx}(dmrsAntIndices{gnbIdx}) + dmrsSymbols{gnbIdx}(:,p)*precodingWeightsIdeal{gnbIdx}(p,:);
        end
      end

      %% OFDM modulate and channel propagation
      % OFDM Modulation - OFDM-modulate the resource grid
      txWaveform = cell(1,num_gNB);
      txWaveformInt = cell(1,num_gNB);
      txWaveformNoInt = cell(1,num_gNB);
      waveformInfo = cell(1,num_gNB);
      for gnbIdx = 1:num_gNB
        [txWaveform{gnbIdx},waveformInfo{gnbIdx}] = nrOFDMModulate(carrier,pdschGrid{gnbIdx});
        [txWaveformInt{gnbIdx},~] = nrOFDMModulate(carrier,pdschGridInt{gnbIdx});
        [txWaveformNoInt{gnbIdx},~] = nrOFDMModulate(carrier,pdschGridNoInt{gnbIdx});
      end

      chInfo = cell(num_gNB,num_UE);
      maxChDelay = zeros(num_gNB,num_UE);
      for gnbIdx = 1:num_gNB
        for ueIdx = 1:num_UE
          % Propagation Channel
          chInfo{gnbIdx,ueIdx} = info(channel{gnbIdx,ueIdx});
          maxChDelay(gnbIdx,ueIdx) = ceil(max(chInfo{gnbIdx,ueIdx}.PathDelays*channel{gnbIdx,ueIdx}.SampleRate)) + chInfo{gnbIdx,ueIdx}.ChannelFilterDelay;
        end
      end
      overallMaxChDelay = max(max(maxChDelay));

      % Pad the input signal with enough zeros to ensure that the generated signal is flushed out of the channel filter.
      for gnbIdx = 1:num_gNB
        txWaveform{gnbIdx} = [txWaveform{gnbIdx}; zeros(overallMaxChDelay,size(txWaveform{gnbIdx},2))];
        txWaveformInt{gnbIdx} = [txWaveformInt{gnbIdx}; zeros(overallMaxChDelay,size(txWaveformInt{gnbIdx},2))];
        txWaveformNoInt{gnbIdx} = [txWaveformNoInt{gnbIdx}; zeros(overallMaxChDelay,size(txWaveformNoInt{gnbIdx},2))];
      end

      % Send the signal through all channels (gNB gnbIdx to UE ueIdx)
      rxWaveform = cell(length(gNB_numTxAnts),length(UE_numRxAnts));
      rxWaveformInt = cell(length(gNB_numTxAnts),length(UE_numRxAnts));
      rxWaveformNoInt = cell(length(gNB_numTxAnts),length(UE_numRxAnts));
      pathGains = cell(length(gNB_numTxAnts),length(UE_numRxAnts));
      pathGainsInt = cell(length(gNB_numTxAnts),length(UE_numRxAnts));
      pathGainsNoInt = cell(length(gNB_numTxAnts),length(UE_numRxAnts));
      sampleTimes = cell(length(gNB_numTxAnts),length(UE_numRxAnts));
      sampleTimesInt = cell(length(gNB_numTxAnts),length(UE_numRxAnts));
      sampleTimesNoInt = cell(length(gNB_numTxAnts),length(UE_numRxAnts));
      for gnbIdx = 1:num_gNB
        for ueIdx = 1:num_UE
          % Propagation Channel
          [rxWaveform{gnbIdx,ueIdx},pathGains{gnbIdx,ueIdx},sampleTimes{gnbIdx,ueIdx}] = channel{gnbIdx,ueIdx}(txWaveform{gnbIdx});
          [rxWaveformInt{gnbIdx,ueIdx},pathGainsInt{gnbIdx,ueIdx},sampleTimesInt{gnbIdx,ueIdx}] = channel{gnbIdx,ueIdx}(txWaveformInt{gnbIdx});
          [rxWaveformNoInt{gnbIdx,ueIdx},pathGainsNoInt{gnbIdx,ueIdx},sampleTimesNoInt{gnbIdx,ueIdx}] = channel{gnbIdx,ueIdx}(txWaveformNoInt{gnbIdx});
        end
      end

      % IAC / Interference / No interference ueRxWaveform
      ueRxWaveform = cell(1,num_UE);
      ueRxWaveformInt = cell(1,num_UE);
      ueRxWaveformNoInt = cell(1,num_UE);
      ueNoiseVec = cell(1,num_UE);
      for ueIdx = 1:num_UE
        % IAC waveform
        ueRxWaveform{ueIdx} = zeros(size(rxWaveform{ueIdx,ueIdx}));
        % Int waveform
        ueRxWaveformInt{ueIdx} = zeros(size(rxWaveformInt{ueIdx,ueIdx}));
        % No Int waveform
        ueRxWaveformNoInt{ueIdx} = rxWaveformNoInt{ueIdx,ueIdx};
        for gnbIdx = 1:num_gNB
          % Combined interfering waveforms for IAC / Int waveforms
          ueRxWaveform{ueIdx} = ueRxWaveform{ueIdx} + rxWaveform{gnbIdx,ueIdx};
          ueRxWaveformInt{ueIdx} = ueRxWaveformInt{ueIdx} + rxWaveformInt{gnbIdx,ueIdx};
        end
        ueNoiseVec{ueIdx} = generateAWGN(SNRdB,UE_numRxAnts(ueIdx),waveformInfo{ueIdx}.Nfft,size(rxWaveform{ueIdx,ueIdx}));
        if noiseless
          ueNoiseVec{ueIdx} = complex(zeros(size(ueNoiseVec{ueIdx})));
        end
        ueRxWaveform{ueIdx} = ueRxWaveform{ueIdx} + ueNoiseVec{ueIdx};
        ueRxWaveformInt{ueIdx} = ueRxWaveformInt{ueIdx} + ueNoiseVec{ueIdx};
        ueRxWaveformNoInt{ueIdx} = ueRxWaveformNoInt{ueIdx} + ueNoiseVec{ueIdx};
      end
      
      %{
      if plotFigures
        pout = fftshift(pwelch(ueRxWaveform{1}));
        figure;plot(10*log10(pout));
        pout = fftshift(pwelch(ueRxWaveform{2}));
        figure;plot(10*log10(pout));
      end
      %}

      %% Decode at UE1, UE2, and UE3
      % Timing Synchronization - Perform timing estimation and synchronization
      pathFilters = cell(num_gNB,num_UE);
      for gnbIdx = 1:num_gNB
        for ueIdx = 1:num_UE
          pathFilters{gnbIdx,ueIdx} = getPathFilters(channel{gnbIdx,ueIdx});
        end
      end

      offset = cell(1,num_UE);
      mag = cell(1,num_UE);
      for ueIdx = 1:num_UE
        [offset{ueIdx},mag{ueIdx}] = nrPerfectTimingEstimate(pathGains{ueIdx,ueIdx},pathFilters{ueIdx,ueIdx});
        % IAC
        ueRxWaveform{ueIdx} = ueRxWaveform{ueIdx}(1+offset{ueIdx}:end,:);
        % Interference
        ueRxWaveformInt{ueIdx} = ueRxWaveformInt{ueIdx}(1+offset{ueIdx}:end,:);
        % No interference
        ueRxWaveformNoInt{ueIdx} = ueRxWaveformNoInt{ueIdx}(1+offset{ueIdx}:end,:);
      end

      % plot waveforms
      %{
      if plotFigures
        figure;
        plot(real(ueRxWaveform{1}),'b');hold on;
        plot(imag(ueRxWaveform{1}),'r');grid on;
        ylim([-0.05 0.05]);title('UE1 RX');legend('real','imag');
        figure;
        subplot(2,1,1);
        plot(real(ueRxWaveform{2}(:,1)),'b');hold on;
        plot(imag(ueRxWaveform{2}(:,1)),'r');grid on;
        ylim([-0.1 0.1]);title('UE2 antenna 1 RX');legend('real','imag');
        subplot(2,1,2);
        plot(real(ueRxWaveform{2}(:,2)),'b');hold on;
        plot(imag(ueRxWaveform{2}(:,2)),'r');grid on;
        ylim([-0.1 0.1]);title('UE2 antenna 2 RX');legend('real','imag');
        figure;
        subplot(3,1,1);
        plot(real(ueRxWaveform{3}(:,1)),'b');hold on;
        plot(imag(ueRxWaveform{3}(:,1)),'r');grid on;
        ylim([-0.15 0.15]);title('UE2 antenna 1 RX');legend('real','imag');
        subplot(3,1,2);
        plot(real(ueRxWaveform{3}(:,2)),'b');hold on;
        plot(imag(ueRxWaveform{3}(:,2)),'r');grid on;
        ylim([-0.15 0.15]);title('UE2 antenna 2 RX');legend('real','imag');
        subplot(3,1,3);
        plot(real(ueRxWaveform{3}(:,3)),'b');hold on;
        plot(imag(ueRxWaveform{3}(:,3)),'r');grid on;
        ylim([-0.15 0.15]);title('UE2 antenna 3 RX');legend('real','imag');
      end
      %}

      % OFDM Demodulation - OFDM-demodulate the synchronized signal
      rxGrid = cell(1,num_UE);
      rxGridInt = cell(1,num_UE);
      rxGridNoInt = cell(1,num_UE);
      for ueIdx = 1:num_UE
        % IAC
        rxGrid{ueIdx} = nrOFDMDemodulate(carrier,ueRxWaveform{ueIdx});
        % Interference
        rxGridInt{ueIdx} = nrOFDMDemodulate(carrier,ueRxWaveformInt{ueIdx});
        % No interference
        rxGridNoInt{ueIdx} = nrOFDMDemodulate(carrier,ueRxWaveformNoInt{ueIdx});
      end

      % Perform perfect channel estimation between transmit and receive antennas.
      estChGridAnts = cell(num_gNB,num_UE);
      estChGridAntsInt = cell(num_gNB,num_UE);
      estChGridAntsNoInt = cell(num_gNB,num_UE);
      noiseGrid = cell(1,num_UE);
      noiseEst = cell(1,num_UE);
      for ueIdx = 1:num_UE
        for gnbIdx = 1:num_gNB
          estChGridAnts{gnbIdx,ueIdx} = nrPerfectChannelEstimate(carrier,pathGains{gnbIdx,ueIdx},pathFilters{gnbIdx,ueIdx},offset{gnbIdx},sampleTimes{gnbIdx,ueIdx});
          estChGridAntsInt{gnbIdx,ueIdx} = nrPerfectChannelEstimate(carrier,pathGainsInt{gnbIdx,ueIdx},pathFilters{gnbIdx,ueIdx},offset{gnbIdx},sampleTimesInt{gnbIdx,ueIdx});
          estChGridAntsNoInt{gnbIdx,ueIdx} = nrPerfectChannelEstimate(carrier,pathGainsNoInt{gnbIdx,ueIdx},pathFilters{gnbIdx,ueIdx},offset{gnbIdx},sampleTimesNoInt{gnbIdx,ueIdx});
        end
        % Get perfect noise estimate (from noise realization)
        noiseGrid{ueIdx} = nrOFDMDemodulate(carrier,ueNoiseVec{ueIdx}(1+offset{ueIdx}:end ,:));
        noiseEst{ueIdx} = var(noiseGrid{ueIdx}(:));
      end

      % Apply precoding to estChGridAnts. The resulting estimate is for the channel estimate between layers and receive antennas.
      estChGridLayers{1,1} = precodeChannelEstimatePerSc(estChGridAnts{1,1},precodingWeights{1}.');
      estChGridLayers{2,1} = estChGridAnts{2,1};
      estChGridLayers{2,2} = estChGridAnts{2,2};

      % Get precoding matrices for next slot
      % obtain precoding weights for gNB1 (align x_1 with x_2 to allow decoding of x_3 at UE2)
      precodingWeights{1} = getPrecodingMatrixIA(pdsch{1}.PRBSet,estChannelGrid{1,2},estChannelGrid{2,2});

      % obtain precoding weights for gNB2
%       precodingWeights{2} = getPrecodingMatrixBeamNull(pdsch{2}.PRBSet,estChGridAnts{2,1});


      %% IAC Equalization
      [pdschRx{1},pdschHest{1,1}] = nrExtractResources(pdschIndices{1},rxGrid{1},estChGridLayers{1,1});
      [~,pdschHest{2,1}] = nrExtractResources(pdschIndices{1},rxGrid{1},estChGridLayers{2,1});
      [pdschRx{2},pdschHest{2,2}] = nrExtractResources(pdschIndices{2},rxGrid{2},estChGridLayers{2,2});

      % solve for (gamma*x1 + x2) and x3 estimate
      pdschEqTemp{2} = complex(zeros(size(pdschRx{2})));
      for i = 1:length(pdschRx{2})
        H = reshape(pdschHest{2,2}(i,:,:),size(pdschHest{2,2},2),size(pdschHest{2,2},3));
        pdschEqTemp{2}(i,:) = pinv(H) * pdschRx{2}(i,:).';
      end
      pdschEqTemp{2} = pdschEqTemp{2}(:,2); % drop gamma*x1 + x2 estimate

      % solve for x1 and x2 estimate
      pdschEqTemp{1} = complex(zeros(size(pdschRx{1})));
      for i = 1:length(pdschRx{1})
        H = [pdschHest{1,1}(i,:).', pdschHest{2,1}(i,:,1).'];
        %     pdschEq_1_2(i,:) = pinv(H) * (pdschRx_1(i,:).' - pdschHest_12(i,:,2).' * pdschEq_3(i));
        pdschEqTemp{1}(i,:) = pinv(H) * (pdschRx{1}(i,:).' - pdschHest{2,1}(i,:,2).' * pdschEqTemp{2}(i));
      end

      [~,csi{1}] = nrEqualizeMMSE(pdschRx{1},pdschHest{1,1},noiseEst{1});
      [~,csi{2}] = nrEqualizeMMSE(pdschRx{2},pdschHest{2,2},noiseEst{2});

      pdschEq{1} = pdschEqTemp{1}(:,1);
      pdschEq{2} = [pdschEqTemp{1}(:,2),pdschEqTemp{2}];


      %% MMSE Equalization (Int / No Int)
      % equalize UE 1,2
      estChGridLayersInt = cell(num_UE,num_UE);
      estChGridLayersNoInt = cell(num_UE,num_UE);
      pdschHestInt = cell(num_UE,num_UE);
      pdschHestNoInt = cell(num_UE,num_UE);
      pdschRxInt = cell(1,num_UE);
      pdschRxNoInt = cell(1,num_UE);
      pdschEqInt = cell(1,num_UE);
      pdschEqNoInt = cell(1,num_UE);
      csiInt = cell(1,num_UE);
      csiNoInt = cell(1,num_UE);
      for ueIdx = 1:num_UE
        % Int
        estChGridLayersInt{ueIdx,ueIdx} = precodeChannelEstimate(estChGridAntsInt{ueIdx,ueIdx},precodingWeightsIdeal{ueIdx}.');
        [pdschRxInt{ueIdx},pdschHestInt{ueIdx,ueIdx}] = nrExtractResources(pdschIndices{ueIdx},rxGridInt{ueIdx},estChGridLayersInt{ueIdx,ueIdx});
        [pdschEqInt{ueIdx},csiInt{ueIdx}] = nrEqualizeMMSE(pdschRxInt{ueIdx},pdschHestInt{ueIdx,ueIdx},noiseEst{ueIdx});
        % No Int
        estChGridLayersNoInt{ueIdx,ueIdx} = precodeChannelEstimate(estChGridAntsNoInt{ueIdx,ueIdx},precodingWeightsIdeal{ueIdx}.');
        [pdschRxNoInt{ueIdx},pdschHestNoInt{ueIdx,ueIdx}] = nrExtractResources(pdschIndices{ueIdx},rxGridNoInt{ueIdx},estChGridLayersNoInt{ueIdx,ueIdx});
        [pdschEqNoInt{ueIdx},csiNoInt{ueIdx}] = nrEqualizeMMSE(pdschRxNoInt{ueIdx},pdschHestNoInt{ueIdx,ueIdx},noiseEst{ueIdx});
      end

      % PDSCH Decoding
      dlschLLRs = cell(1,num_UE);
      dlschLLRsInt = cell(1,num_UE);
      dlschLLRsNoInt = cell(1,num_UE);
      rxSymbols = cell(1,num_UE);
      rxSymbolsInt = cell(1,num_UE);
      rxSymbolsNoInt = cell(1,num_UE);
      for ueIdx = 1:num_UE
        if plotFigures
          figure;plot(real(pdschEq{ueIdx}),imag(pdschEq{ueIdx}),'b.');grid on;
          xlabel('In-Phase');ylabel('Quadrature');title('UE' + string(ueIdx) + ' IAC RX IQ');axis([-1.5 1.5 -1.5 1.5]);axis square;
          figure;plot(real(pdschEqInt{ueIdx}),imag(pdschEqInt{ueIdx}),'b.');grid on;
          xlabel('In-Phase');ylabel('Quadrature');title('UE' + string(ueIdx) + ' Int RX IQ');axis([-1.5 1.5 -1.5 1.5]);axis square;
          figure;plot(real(pdschEqNoInt{ueIdx}),imag(pdschEqNoInt{ueIdx}),'b.');grid on;
          xlabel('In-Phase');ylabel('Quadrature');title('UE' + string(ueIdx) + ' No Int RX IQ');axis([-1.5 1.5 -1.5 1.5]);axis square;
        end

        %% IAC Decode
        [dlschLLRs{ueIdx},rxSymbols{ueIdx}] = nrPDSCHDecode(carrier,pdsch{ueIdx},pdschEq{ueIdx},noiseEst{ueIdx});

        % Scale LLRs by CSI
        csi{ueIdx} = nrLayerDemap(csi{ueIdx}); % CSI layer demapping
        for cwIdx = 1:pdsch{ueIdx}.NumCodewords
          Qm = length(dlschLLRs{ueIdx}{cwIdx})/length(rxSymbols{ueIdx}{cwIdx}); % Bits per symbol
          csi{ueIdx}{cwIdx} = repmat(csi{ueIdx}{cwIdx}.',Qm,1); % Expand by each bit per symbol
          dlschLLRs{ueIdx}{cwIdx} = dlschLLRs{ueIdx}{cwIdx} .* csi{ueIdx}{cwIdx}(:); % Scale
        end

        decodeDLSCH{ueIdx}.TransportBlockLength = trBlkSizes{ueIdx};
        [decbits{nSlot+1,ueIdx},blkerr{nSlot+1,ueIdx}] = decodeDLSCH{ueIdx}(dlschLLRs{ueIdx},pdsch{ueIdx}.Modulation,pdsch{ueIdx}.NumLayers, ...
          harqEntity{ueIdx}.RedundancyVersion,harqEntity{ueIdx}.HARQProcessID);


        %% Int Decode
        [dlschLLRsInt{ueIdx},rxSymbolsInt{ueIdx}] = nrPDSCHDecode(carrier,pdsch{ueIdx},pdschEqInt{ueIdx},noiseEst{ueIdx});

        % Scale LLRs by CSI
        csiInt{ueIdx} = nrLayerDemap(csiInt{ueIdx}); % CSI layer demapping
        for cwIdx = 1:pdsch{ueIdx}.NumCodewords
          Qm = length(dlschLLRsInt{ueIdx}{cwIdx})/length(rxSymbolsInt{ueIdx}{cwIdx}); % Bits per symbol
          csiInt{ueIdx}{cwIdx} = repmat(csiInt{ueIdx}{cwIdx}.',Qm,1); % Expand by each bit per symbol
          dlschLLRsInt{ueIdx}{cwIdx} = dlschLLRsInt{ueIdx}{cwIdx} .* csiInt{ueIdx}{cwIdx}(:); % Scale
        end

        decodeDLSCHInt{ueIdx}.TransportBlockLength = trBlkSizes{ueIdx};
        [decbitsInt{nSlot+1,ueIdx},blkerrInt{nSlot+1,ueIdx}] = decodeDLSCHInt{ueIdx}(dlschLLRsInt{ueIdx},pdsch{ueIdx}.Modulation,pdsch{ueIdx}.NumLayers, ...
          harqEntityInt{ueIdx}.RedundancyVersion,harqEntityInt{ueIdx}.HARQProcessID);

        
        %% No Int Decode
        [dlschLLRsNoInt{ueIdx},rxSymbolsNoInt{ueIdx}] = nrPDSCHDecode(carrier,pdsch{ueIdx},pdschEqNoInt{ueIdx},noiseEst{ueIdx});

        % Scale LLRs by CSI
        csiNoInt{ueIdx} = nrLayerDemap(csiNoInt{ueIdx}); % CSI layer demapping
        for cwIdx = 1:pdsch{ueIdx}.NumCodewords
          Qm = length(dlschLLRsNoInt{ueIdx}{cwIdx})/length(rxSymbolsNoInt{ueIdx}{cwIdx}); % Bits per symbol
          csiNoInt{ueIdx}{cwIdx} = repmat(csiNoInt{ueIdx}{cwIdx}.',Qm,1); % Expand by each bit per symbol
          dlschLLRsNoInt{ueIdx}{cwIdx} = dlschLLRsNoInt{ueIdx}{cwIdx} .* csiNoInt{ueIdx}{cwIdx}(:); % Scale
        end

        decodeDLSCHNoInt{ueIdx}.TransportBlockLength = trBlkSizes{ueIdx};
        [decbitsNoInt{nSlot+1,ueIdx},blkerrNoInt{nSlot+1,ueIdx}] = decodeDLSCHNoInt{ueIdx}(dlschLLRsNoInt{ueIdx},pdsch{ueIdx}.Modulation,pdsch{ueIdx}.NumLayers, ...
          harqEntityNoInt{ueIdx}.RedundancyVersion,harqEntityNoInt{ueIdx}.HARQProcessID);


        %% Decode TX bits
        harqEntityCpy = harqEntity{ueIdx};
        [txLlrs,~] = nrPDSCHDecode(carrier,pdsch{ueIdx},pdschSymbolsRaw{ueIdx},0);
        [txBits,~] = decodeDLSCH{ueIdx}(txLlrs,pdsch{ueIdx}.Modulation,pdsch{ueIdx}.NumLayers, ...
          harqEntityCpy.RedundancyVersion,harqEntityCpy.HARQProcessID);
        harqEntityCpyInt = harqEntityInt{ueIdx};
        [txLlrsInt,~] = nrPDSCHDecode(carrier,pdsch{ueIdx},pdschSymbolsRawInt{ueIdx},0);
        [txBitsInt,~] = decodeDLSCHInt{ueIdx}(txLlrsInt,pdsch{ueIdx}.Modulation,pdsch{ueIdx}.NumLayers, ...
          harqEntityCpyInt.RedundancyVersion,harqEntityCpyInt.HARQProcessID);
        harqEntityCpyNoInt = harqEntity{ueIdx};
        [txLlrsNoInt,~] = nrPDSCHDecode(carrier,pdsch{ueIdx},pdschSymbolsRawNoInt{ueIdx},0);
        [txBitsNoInt,~] = decodeDLSCHNoInt{ueIdx}(txLlrsNoInt,pdsch{ueIdx}.Modulation,pdsch{ueIdx}.NumLayers, ...
          harqEntityCpyNoInt.RedundancyVersion,harqEntityCpyNoInt.HARQProcessID);

        berIac(ueIdx,nSlot+1)   = berIac(ueIdx,nSlot+1)   + sum(~(decbits{nSlot+1,ueIdx}      == txBits)) / length(txBits); % BER
        berInt(ueIdx,nSlot+1)   = berInt(ueIdx,nSlot+1)   + sum(~(decbitsInt{nSlot+1,ueIdx}   == txBitsInt)) / length(txBitsInt); % BER
        berNoInt(ueIdx,nSlot+1) = berNoInt(ueIdx,nSlot+1) + sum(~(decbitsNoInt{nSlot+1,ueIdx} == txBitsNoInt)) / length(txBitsNoInt); % BER
      end
    end
  end
  BER_IAC(snrIdx,:,:) = berIac / numMC;
  BER_INT(snrIdx,:,:) = berInt / numMC;
  BER_NO_INT(snrIdx,:,:) = berNoInt / numMC;
  mc_loop_time = toc(simTimer);
  fprintf('SNR index %i out of %i completed at %0.2f seconds\n',snrIdx,length(snrRange),mc_loop_time);
end

figure;
semilogy(snrRange, BER_IAC(:,1),'b.-','linewidth',1.2);hold on;semilogy(snrRange, BER_NO_INT(:,1),'r.-','linewidth',1.2);semilogy(snrRange, BER_INT(:,1),'k.-','linewidth',1.2);grid on;
title('Average BER');xlabel('SNR (dB)');ylabel('Average BER');axis square;
legend('IAC','No Interference','Interference');title('UE1 performance');
ylim([1e-3 1]);

figure;
semilogy(snrRange, BER_IAC(:,2),'b.-','linewidth',1.2);hold on;semilogy(snrRange, BER_NO_INT(:,2),'r.-','linewidth',1.2);semilogy(snrRange, BER_INT(:,2),'k.-','linewidth',1.2);grid on;
title('Average BER');xlabel('SNR (dB)');ylabel('Average BER');axis square;
legend('IAC','No Interference','Interference');title('UE2 performance');
ylim([1e-3 1]);


%% Local Functions
%% generateAWGN
function noise = generateAWGN(SNRdB,nRxAnts,Nfft,sizeRxWaveform)
  SNR = 10^(SNRdB/10); % Calculate linear noise gain
  N0 = 1/sqrt(2.0*nRxAnts*double(Nfft)*SNR);
  noise = N0*complex(randn(sizeRxWaveform),randn(sizeRxWaveform));
end


%% getPrecodingMatrixSvd
function svd_wtx = getPrecodingMatrixSvd(PRBSet,NLayers,hestGrid)
  % Calculate precoding matrix given an allocation and a channel estimate

  % Allocated subcarrier indices
  allocSc = (1:12)' + 12*PRBSet(:).';
  allocSc = allocSc(:);

  % Average channel estimate
  [~,~,R,P] = size(hestGrid);
  estAllocGrid = hestGrid(allocSc,:,:,:); % extact allocated symbols (keeps all OFDM symbols and TX/RX antennas)
  Hrs = reshape(estAllocGrid,[],R,P); % stack all OFDM symbols in the first dimension
  Havg = mean(Hrs); % averaging all subcarriers accross all OFDM symbols
  Hest = permute(Havg,[2 3 1]); % drop first dimension

  % SVD decomposition
  [~,~,V] = svd(Hest); % using SVD to compute the precoding matrix (V :: N x N complex unitary matrix where N = Num TX antennas)

  svd_wtx = V(:,1:NLayers).';
  svd_wtx = svd_wtx/sqrt(NLayers); % Normalize by NLayers
end


%% getPrecodingMatrixIA
function ia_weights = getPrecodingMatrixIA(PRBSet,hestGrid_12,hestGrid_22)
  % Calculate precoding matrix given an allocation and a channel estimate

  % Allocated subcarrier indices
  allocSc = (1:12)' + 12*PRBSet(:).';
  allocSc = allocSc(:);

  % precode per subcarrier
  [~,~,R_12,P_12] = size(hestGrid_12);
  estAllocGrid_12 = hestGrid_12(allocSc,:,:,:); % extact allocated symbols (keeps all OFDM symbols and TX/RX antennas)
  [~,~,R_22,P_22] = size(hestGrid_22);
  estAllocGrid_22 = hestGrid_22(allocSc,:,:,:); % extact allocated symbols (keeps all OFDM symbols and TX/RX antennas)

  Havg_12 = mean(estAllocGrid_12,2); % average over OFDM symbols
  Havg_22 = mean(estAllocGrid_22,2); % average over OFDM symbols
  Hest_12 = reshape(Havg_12,[],R_12,P_12); % drop second dimension
  Hest_22 = reshape(Havg_22,[],R_22,P_22); % drop second dimension

  % do the computation per-sc loop and then do it flattened (compare results, if they match... drop the loop)
  ia_weights = complex(zeros(size(Hest_12,1), size(Hest_12,2)));
  for scIdx = 1:size(Hest_12,1)
    H = reshape(Hest_12(scIdx,:,:),R_22,P_22);
    ia_weights(scIdx,:) = pinv(H) * Hest_22(scIdx,:,1).';
  end
end


%% getInitialChannelEstimate
function estChannelGrid = getInitialChannelEstimate(channel,carrier)
  % Obtain an initial channel estimate for calculating the precoding matrix.
  % This function assumes a perfect channel estimate

  % Clone of the channel
  chClone = channel.clone();
  chClone.release();

  % No filtering needed to get channel path gains
  chClone.ChannelFiltering = false;

  % Get channel path gains
  [pathGains,sampleTimes] = chClone();

  % Perfect timing synchronization
  pathFilters = getPathFilters(chClone);
  offset = nrPerfectTimingEstimate(pathGains,pathFilters);

  % Perfect channel estimate
  estChannelGrid = nrPerfectChannelEstimate(carrier,pathGains,pathFilters,offset,sampleTimes);
end


%% precodeChannelEstimate
function estChannelGrid = precodeChannelEstimate(estChannelGrid,W)
  K = size(estChannelGrid,1); % num subcarriers
  L = size(estChannelGrid,2); % num ofdm symbols
  R = size(estChannelGrid,3); % number of RX antennas
  %   T = size(W,2);              % number of TX antennas
  estChannelGrid = reshape(estChannelGrid,K*L*R,[]);
  estChannelGrid = estChannelGrid*W;
  estChannelGrid = reshape(estChannelGrid,K,L,R,[]);
end


%% precodeChannelEstimatePerSc
function pcEstChannelGrid = precodeChannelEstimatePerSc(estChannelGrid,W)
  K = size(estChannelGrid,1); % num subcarriers
  L = size(estChannelGrid,2); % num ofdm symbols
  R = size(estChannelGrid,3); % number of RX antennas
  T = size(estChannelGrid,4); % number of TX antennas

  % gNB 2 beam nulling (loop through symbols)
  pcEstChannelGrid = zeros(K,L,R);
  for symIdx = 1:L
    %     symEstChannelGrid = estChannelGrid(:,symIdx,:,:);
    %     pcEstChannelGrid(:,symIdx,:,:) = symEstChannelGrid .* W;
    symEstChannelGrid = reshape(estChannelGrid(:,symIdx,:,:),K,R,T,[]);
    for scIdx = 1:K
      scChannelEst = reshape(symEstChannelGrid(scIdx,:,:),R,T,[]);
      pcEstChannelGrid(scIdx,symIdx,:) = scChannelEst * W(:,scIdx);
    end
  end
end

