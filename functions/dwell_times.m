function [ times ] = dwell_times(seq,s)
%DWELL_TIMES Given a sequence of states (seq),
% compute the dwell times state s in the sequence.
    imax = length(seq);
    times = [];
    previous_state = seq(1);
    counter = 0;
    for idx = 2:imax
     if(previous_state == seq(idx))
     counter = counter + 1;
     else
     if (previous_state == s)
     times = [times ; counter];
     end
     counter = 1;
     previous_state = seq(idx);
     end
    end

end