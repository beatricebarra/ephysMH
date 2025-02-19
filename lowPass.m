function signals = lowPass(signals, fS, fCut, order)
	% This method returns a low-passed version of the input signals.
	%
	% Inputs: - signals: N*nChan double-array: nChan signals.
	%		  - fS: double: sampling frequency.
	%		  - fCut: double: cutting frequency.
	%		  - order: integer: filter order.
	%
	% Outputs - signals: N*nChan double-array: the low-passed signals.
	%
	% Actions: None.
	%
	% Description: Void.

	% Parse inputs.
	if (nargin == 3)
		order = 3;
	end

	fC = 2*fCut/fS;
	[N, nChan] = size(signals);

	% Low-pass filtering.
	d = fdesign.lowpass('N,Fc', order, fC);
	h = design(d, 'butter');
	for(i=1:nChan)
	    signals(:,i) = filtfilt(h.sosMatrix, h.ScaleValues, signals(:,i)');
	end

end