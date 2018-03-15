
connectedRadios = findsdru;
if strncmp(connectedRadios(1).Status, 'Success', 7)
  radioFound = true;
  platform = connectedRadios(1).Platform;
  switch connectedRadios(1).Platform
    case {'B200','B210'}
      address = connectedRadios(1).SerialNum;
    case {'N200/N210/USRP2','X300','X310'}
      address = connectedRadios(1).IPAddress;
  end
else
  radioFound = false;
  address = '192.168.10.2';
  platform = 'N200/N210/USRP2';
end

fmRxParams = getParamsSdruFMExamples(platform)

% Set up radio object to use the found radio
switch platform
  case {'B200','B210'}
    radio = comm.SDRuReceiver(...
      'Platform', platform, ...
      'SerialNum', address, ...
      'MasterClockRate', fmRxParams.RadioMasterClockRate);
  case {'X300','X310'}
    radio = comm.SDRuReceiver(...
      'Platform', platform, ...
      'IPAddress', address, ...
      'MasterClockRate', fmRxParams.RadioMasterClockRate);
  case {'N200/N210/USRP2'}
    radio = comm.SDRuReceiver(...
      'Platform', platform, ...
      'IPAddress', address);
end
radio.CenterFrequency  = 102.5e6;
radio.Gain = fmRxParams.RadioGain;
radio.DecimationFactor = fmRxParams.RadioDecimationFactor;
radio.SamplesPerFrame = fmRxParams.RadioFrameLength;
radio.OutputDataType = 'single'
hwInfo = info(radio)
fmBroadcastDemod = comm.FMBroadcastDemodulator(...
    'SampleRate', fmRxParams.RadioSampleRate, ...
    'FrequencyDeviation', fmRxParams.FrequencyDeviation, ...
    'FilterTimeConstant', fmRxParams.FilterTimeConstant, ...
    'AudioSampleRate', fmRxParams.AudioSampleRate, ...
    'PlaySound', true, ...
    'BufferSize', fmRxParams.BufferSize, ...
    'Stereo', true);
if radioFound
  % Loop until the example reaches the target stop time, which is 10
  % seconds.
  timeCounter = 0;
  while timeCounter < fmRxParams.StopTime
    [x, len] = step(radio);
    if len > 0
      % FM demodulation
      step(fmBroadcastDemod, x);
      % Update counter
      timeCounter = timeCounter + fmRxParams.AudioFrameTime;
    end
  end
else
  warning(message('sdru:sysobjdemos:MainLoop'))
end

release(fmBroadcastDemod)
release(radio)