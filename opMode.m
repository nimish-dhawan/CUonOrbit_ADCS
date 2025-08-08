function MaskInitialization(maskInitContext)

    % Get the current value of the mask parameter 'opMode'
    opModeValue = maskInitContext.MaskObject.get('opMode');

    % Create a Simulink.Parameter object to be used as a variant control
    OpMode = Simulink.Parameter;
    OpMode.Value = str2double(opModeValue);  % Convert to number if stored as string
    OpMode.CoderInfo.StorageClass = 'Auto';

    % Push the variable to base workspace (for use in variant control expressions)
    assignin('base', 'opMode', OpMode);

end
