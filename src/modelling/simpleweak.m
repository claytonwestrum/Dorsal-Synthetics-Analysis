function yfit = simpleweak(x, params, varargin)
%simplebinding in the weak promoter limit.

noOff = false;
fraction = false;
dimer = false;


%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end


if isstruct(x)
    data = x;
    x = data.X(:, 1);
end

X = [];
n = 1;
X(1, :) = [x(1), 1];
for k = 2:length(x)
    if x(k) < x(k-1)
        n = n + 1;
    end
    X(k, 1) = x(k);
    X(k, 2) = n;
end

dsid = X(:,2);     % unpack dataset id from X
params = params(:)'; %need a row vec

if ~noOff
    offset = params(max(dsid)+3);
else
    offset = 0;
end

if ~fraction
    amp = params(max(dsid)+2);
else
    amp = 1;
    offset = 0;
end

if expmnt=="affinities"
    omegaDP = params(1);
    KD = params(2:max(dsid)+1)';
    KD = KD(dsid);
elseif expmnt == "phases"
    KD = params(1);
    omegaDP = params(2:max(dsid)+1)';
    omegaDP = omegaDP(dsid);
elseif expmnt == "phaff"
    nph = 3;
    naff = 7;
    omegaDP = [params(1)*ones(1, naff), params(2:nph)]';
    omegaDP = omegaDP(dsid);
    KD = [params(1:naff), params(1)*ones(1, nph-1)]';
    KD = KD(dsid);
    if ~fraction
        amp = params(nph+naff+1);
        offset = 0;
    end
end


if dimer
    k = params(end); %dimerization kd
    x = (1/2).*(k+2.*x+(-1).*k.^(1/2).*(k+4.*x).^(1/2)); %concentration of dorsal dimers
end

sb = @(x, p) p{1}.*(((x./p{2}).*p{3})./(1+ (x./p{2})+ ((x./p{2}).*p{3}))) + p{4};
%sb inputs are 1. dorsal 2. {amp, KD, omegaDP, offset}. Any of the elements
%of 2 are allowed to be column vectors.
yfit = sb(x, {amp; KD; omegaDP; offset});

end