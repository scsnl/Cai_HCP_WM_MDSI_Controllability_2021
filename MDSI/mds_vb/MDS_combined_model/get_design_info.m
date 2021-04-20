function [Vm, names] = get_design_info (T,TR,spm_mat,sess_number)
load(spm_mat);
num_events = length(SPM.Sess(sess_number).U);
names = {};
onsets = {};
durations = {};
Vm = zeros(T,num_events);

for i = 1:num_events
  names{i} = SPM.Sess(sess_number).U(i).name{1};
  onsets{i} = SPM.Sess(sess_number).U(i).ons;
  durations{i} = SPM.Sess(sess_number).U(i).dur;
  
  for t = 1:length(onsets{i})
    ns = floor(onsets{i}(t)/TR) + 1; %Start of event
    Vm(ns,i) = 1;
    ne = ns + floor(durations{i}(t)/TR)-1; % End of the event;
    Vm(ns:ne,i) = 1;
  end
end

if ~isempty(find(sum(Vm,2) > 1, 1))
  error('Task Design is not orthogonal');
elseif sum(sum(Vm,2)) < size(Vm,1)
  names{num_events+1} = 'rest';
  Vm(:,num_events+1) = ones(size(Vm,1),1) - sum(Vm,2);
end
 
clear SPM;

end
