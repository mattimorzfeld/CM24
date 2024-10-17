function [start,stop]=FindEventsActivity(q,ST,RT)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% t : Time 
% q1 : quantity 1 (PSV index)
% q2 : quantity 2 (dipole or PSV index)
% ST: Start of event threshold (event starts when q EXCEEDS ST)
% RT: Recovery threshold (event stops when q IS BELOW RT)
%
% Output
% start: Time index at start of each event
% stop : Time index at end of each event
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start=[];
stop=[];

% flag = 1 during event, flag =0 otherwise
if abs(q(1))>ST
    flag=1;
    start=[start;1];
else
    flag=0;
end

for ii=2:length(q)
    if flag==0 && q(ii)>ST
        start=[start;ii];
        flag=1;
    elseif flag==1 && q(ii)<RT
        stop=[stop;ii];
        flag=0;
    end
end

if length(start)~=length(stop)
    start=start(1:end-1);
end

        
        