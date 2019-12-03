function options = spgSetParms(varargin)
%SPGSETPARMS  Set options for SPGL1
%
%   options = spgSetParms('param1',val1,'param2',val2,...) creates an
%   options structure in which the named parameters have the specified
%   values.  Unspecified parameters are empty and their default
%   values are used.
%   
%   spgSetParms with no input arguments and no output arguments
%   displays all parameter names and their possible values.
%
%   options = spgSetParms (with no input arguments) creates an options
%   structure where all the fields are empty.
%
%   spgSetParms.m
%   $Id: spgSetParms.m 1407 2009-06-30 20:00:54Z ewout78 $



parametersClassic = {...
    {'fid',         'positive integer',          '1'    }, ...
    {'verbosity',   'integer: 1, 2, or 3',       '2'    }, ...
    {'history',     'false=no, true=yes',        'false'}, ...
    {'iterations',  'positive integer',          'NaN'  }, ...
    {'nPrevVals',   'positive integer',          '3'    }, ...
    {'bpTol',       'positive scalar',           'NaN'  }, ...
    {'lsTol',       'positive scalar',           '1e-06'}, ...
    {'optTol',      'positive scalar',           '1e-04'}, ...
    {'decTol',      'positive scalar',           '1e-04'}, ...
    {'projTol',     'positive scalar',           'NaN'  }, ...
    {'relgapMinF',  'positive scalar',           '1'    }, ...
    {'relgapMinR',  'positive scalar',           '1'    }, ...
    {'rootfindMode','0 = primal, 1 = dual',      '0'    }, ...
    {'rootfindTol', 'scalar',                    '0.5'  }, ...
    {'stepMin',     'positive scalar',           '1e-16'}, ...
    {'stepMax',     'positive scalar',           '1e+05'}, ...
    {'iscomplex',   '0=no, 1=yes, NaN=auto',     'NaN'  }, ...
    {'maxMatvec',   'positive integer',          'Inf'  }, ...
    {'maxRuntime',  'runtime in seconds',        'Inf'  }, ...
    {'mu',          'nonnegative scalar',        '0'    }, ...
    {'weights',     'vector or scalar',          '1'    }, ...
    {'project',     'projection function',       '@NormL1_project'}, ...
    {'primal_norm', 'primal norm eval fun',      '@NormL1_primal'}, ...
    {'dual_norm',   'dual norm eval fun',        '@NormL1_dual'}, ...
    {'hybridMode',  'false=no, true=yes',        'false'}, ...  % Hybrid
    {'lbfgsHist',   'positive integer',          '8'    }, ...
   };


parametersHybrid = ...
   {{'fid',         'positive integer',          '1'    }, ...
    {'verbosity',   'integer: 1, 2, or 3',       '2'    }, ...
    {'history',     'false=no, true=yes',        'false'}, ...
    {'iterations',  'positive integer',          'NaN'  }, ...
    {'nPrevVals',   'positive integer',          '3'    }, ...
    {'bpTol',       'positive scalar',           '0'    }, ...
    {'lsTol',       'positive scalar',           '0'    }, ...
    {'optTol',      'positive scalar',           '1e-04'}, ...
    {'decTol',      'positive scalar',           '1e-04'}, ...
    {'projTol',     'positive scalar',           'NaN'  }, ...
    {'relgapMinF',  'positive scalar',           '1'    }, ...
    {'relgapMinR',  'positive scalar',           '1'    }, ...
    {'rootfindMode','0 = primal, 1 = dual',      '1'    }, ...
    {'rootfindTol', 'scalar',                    '0.9'  }, ...
    {'stepMin',     'positive scalar',           '1e-16'}, ...
    {'stepMax',     'positive scalar',           '1e+05'}, ...
    {'iscomplex',   '0=no, 1=yes, NaN=auto',     'NaN'  }, ...
    {'maxMatvec',   'positive integer',          'Inf'  }, ...
    {'maxRuntime',  'runtime in seconds',        'Inf'  }, ...
    {'mu',          'nonnegative scalar',        '0'    }, ...
    {'weights',     'vector or scalar',          '1'    }, ...
    {'project',     'projection function',       '@NormL1_project'}, ...
    {'primal_norm', 'primal norm eval fun',      '@NormL1_primal'}, ...
    {'dual_norm',   'dual norm eval fun',        '@NormL1_dual'}, ...
    {'hybridMode',  'false=no, true=yes',        'true' }, ...  % Hybrid
    {'lbfgsHist',   'positive integer',          '8'    }, ...
   };


% The first input argument can be a mode
parameters = parametersClassic;
if ((nargin >= 1) && (ischar(varargin{1})))
   nmodes = 1;
   switch varargin{1}
      case {'classic'}
         parameters = parametersClassic;
      case {'hybrid'}
         parameters = parametersHybrid;
      otherwise
         nmodes = 0;
   end
else
   nmodes = 0;
end

% Print out possible values of properties.
if (((nargin-nmodes) == 0) && (nargout == 0))
   fprintf(' Default parameters for spgSetParams.m:\n');
   w1 = 0; w2 = 0; w3 = 0;
   for i = 1:length(parameters)
      w1 = max(w1, length(parameters{i}{1}));
      w2 = max(w2, length(parameters{i}{2}));
      w3 = max(w3, length(parameters{i}{3}));
   end
   for i = 1:length(parameters)
      fprintf('%*s : %*s | %*s\n', w1, parameters{i}{1}, ...
                                  -w2, parameters{i}{2}, ...
                                   w3, parameters{i}{3});
   end
   fprintf('\n');
   return;
end

% Analyze the preset parameters
m = length(parameters);
Names = cell(m,1);
names = cell(m,1);
for i = 1:m
   Names{i} = parameters{i}{1};
   names{i} = lower(Names{i});
end

% Combine all leading options structures o1, o2, ... in l1Set(o1,o2,...).
options = [];
for j = 1:m
   options.(parameters{j}{1}) = eval(parameters{j}{3});
end

i = 1+nmodes;
while i <= nargin
   arg = varargin{i};
   if ischar(arg), break; end
   if ~isempty(arg)                      % [] is a valid options argument
      if ~isa(arg,'struct')
          error(sprintf(['Expected argument %d to be a string parameter name ' ...
               'or an options structure\ncreated with OPTIMSET.'], i));
      end
      for j = 1:m
          if any(strcmp(fieldnames(arg),Names{j}))
             val = arg.(Names{j});
             options.(Names{j}) = val;
          end
      end
   end
   i = i + 1;
end

% A finite state machine to parse name-value pairs.
if rem(nargin-i+1,2) ~= 0
   error('Arguments must occur in name-value pairs.');
end
expectval = 0;                          % start expecting a name, not a value
while i <= nargin
   arg = varargin{i};
   
   if ~expectval
      if ~ischar(arg)
         error(sprintf('Expected argument %d to be a string parameter name.', i));
      end
      
      lowArg = lower(arg);
      j = strmatch(lowArg,names);
      if isempty(j)                       % if no matches
         error(sprintf('Unrecognized parameter name ''%s''.', arg));
      elseif length(j) > 1                % if more than one match
         % Check for any exact matches (in case any names are subsets of others)
         k = strmatch(lowArg,names,'exact');
         if length(k) == 1
            j = k;
         else
            msg = sprintf('Ambiguous parameter name ''%s'' ', arg);
            msg = [msg '(' Names{j(1)}];
            for k = j(2:length(j))'
               msg = [msg ', ' Names{k}];
            end
            msg = sprintf('%s).', msg);
            error(msg);
         end
      end
      expectval = 1;                      % we expect a value next
      
   else
      options.(Names{j}) = arg;
      expectval = 0;
      
   end
   i = i + 1;
end

if expectval
   error(sprintf('Expected value for parameter ''%s''.', arg));
end

