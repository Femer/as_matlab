function [pixTime_us] = myTime2pixhawkTime(myT)
%myTime2pixhawkTime converts my local time (UTC+1) to the pixhawk 
%"absolute" time, given the local time when the pixhaek was switched ON

%now I use the time from my watch at which pixhawk was switched on MINUS
% 2 seconds only to make sure that every timestamp at which I think the
% step command was given, is a little bit BEFORe the real step command. So
% from instead of the original time 10:03:10.00, I use 10:03:12.00
tSwitchOnLocal = '10:03:12.00';

relSec =  datenum(myT, 'HH:MM:SS') .* (24*60*60) - ...
          datenum(tSwitchOnLocal, 'HH:MM:SS') .* (24*60*60);
pixTime_us = relSec * 1e6;
end

