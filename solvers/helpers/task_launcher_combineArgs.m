function varargout = task_launcher_combineArgs(task, varargin )
%task_launcher_combineArgs launches jobs by combining all parameters

if ischar(task)
   type_task = exist(task);
   if type_task~=2 && type_task ~=5
       error('Undefined function %s.m', task);
   end
   task_fn = str2func(task);    
else
   error('First argument should be a string with name of function to be launched'); 
end

% check number of arguments
nargin_fn  = nargin(task_fn);
nargout_fn = nargout(task_fn);

noOutput = false;
if nargout>nargout_fn
    error('Too many output arguments for %s', task);
elseif nargout == 0;
    noOutput = true;
end

if nargin_fn < numel(varargin)
    error('Wrong number of arguments for %s', task);
end

% evaluate number of jobs
n_jobs = prod(cellfun(@numel, varargin));

% combine arguments   
combine = 1:length(varargin{1});

for i = 2 : length(varargin)            
    newcombine = kron(1:length(varargin{i}), ones(1, size(combine, 2)));
    combine = repmat(combine, 1, length(varargin{i}));            
    combine = [combine; newcombine];
end 

% launching jobs
fprintf('Launching %d jobs for %s \n', n_jobs, task);

if ~noOutput
   OUTPUT = cell(1, n_jobs);
   ARGS   = cell(1, n_jobs);
end

for j=1:n_jobs
    cur_arg = cell(1, numel(varargin));
    for a=1:numel(varargin)
        cur_arg{a} = varargin{a}{combine(a,j)};
    end
    % launching task number j
    fprintf('Launching task %d out of %d....\n', j, n_jobs);
    if noOutput
        task_fn(cur_arg{:});
    else
        cur_output = cell(1, nargout);
        [cur_output{1:nargout}] = task_fn(cur_arg{:});
    end
            
    fprintf('Job %d finished successfuly!\n', j);   
    if ~noOutput
       ARGS{j}   = cur_arg;
       OUTPUT{j} = cur_output;
    end
    
end

if ~noOutput
   varargout{1} = OUTPUT;
   varargout{2} = ARGS;   
end

end


