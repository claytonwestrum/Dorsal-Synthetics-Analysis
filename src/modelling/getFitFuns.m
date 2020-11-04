function fun = getFitFuns(expmnt, md, metric, noOff)

if expmnt == "affinities" && contains(md, "simpleweak")
        fun = @(x, p)simpleweak(x, p, 'noOff', noOff, 'fraction', metric=="fraction", 'dimer', contains(md, "dimer"));
elseif expmnt == "phases" && contains(md, "simpleweak")
        fun = @(x, p)simpleweak_phases(x, p, 'noOff', noOff, 'fraction', metric=="fraction", 'dimer', contains(md, "dimer"));
elseif expmnt == "phaff" && contains(md, "simpleweak")
        fun = @(x, p)simpleweak_phaff(x, p, 'noOff', noOff, 'fraction', metric=="fraction", 'dimer', contains(md, "dimer"));
end


end

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
omegaDP = params(1);
KD = params(2:max(dsid)+1)';

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

if dimer
    k = params(end); %dimerization kd 
    x = (1/2).*(k+2.*x+(-1).*k.^(1/2).*(k+4.*x).^(1/2)); %concentration of dorsal dimers
end

%sb inputs are 1. dorsal 2. {amp, KD, omegaDP, offset}. Any of the elements
%of 2 are allowed to be column vectors.
yfit = sb(x, {amp; KD(dsid); omegaDP; offset});

end


function yfit = simpleweak_phases(x, params, varargin)
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
KD = params(1);
omegaDP = params(2:max(dsid)+1)';

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

if dimer
    k = params(end); %dimerization kd 
    x = (1/2).*(k+2.*x+(-1).*k.^(1/2).*(k+4.*x).^(1/2)); %concentration of dorsal dimers
end

%sb inputs are 1. dorsal 2. {amp, KD, omegaDP, offset}. Any of the elements
%of 2 are allowed to be column vectors.
yfit = sb(x, {amp, KD, omegaDP(dsid), offset});
end

function yfit = simpleweak_phaff(x, params, varargin)
%simplebinding in the weak promoter limit.
%note that this could probably be cleaner if dsid was a 2d matrix and
%params were referenced in that matrix like data set 1 has kd1, omega1,
%etc.

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

nph = 3;
naff = 7;
dsid = X(:,2);     % unpack dataset id from X
params = params(:)'; %need a row vec
omegaDP = [params(1)*ones(1, naff), params(2:nph)]';
KD = [params(1:naff), params(1)*ones(1, nph-1)]';

%this is just always nooff
% if ~noOff
%     offset = params(nph+naff+2);
% else
%     offset = 0;
% end

if ~fraction
    amp = params(nph+naff+1);
else
    amp = 1;
end
offset = 0;

if dimer
    k = params(end); %dimerization kd 
    x = (1/2).*(k+2.*x+(-1).*k.^(1/2).*(k+4.*x).^(1/2)); %concentration of dorsal dimers
end

%sb inputs are 1. dorsal 2. {amp, KD, omegaDP, offset}. Any of the elements
%of 2 are allowed to be column vectors.
yfit = sb(x, {amp, KD(dsid), omegaDP(dsid), offset});

end



function m = sb(x, p)

%simple binding model 

%sb inputs are 1. dorsal 2. {amp, KD, omegaDP, offset}. Any of the elements
%of 2 are allowed to be column vectors.

m = p{1}.*(((x./p{2}).*p{3})./(1+ (x./p{2})+ ((x./p{2}).*p{3}))) + p{4};

end