function [pkts,hasData,groupMembers] = baseline(num_users,winner,pkts)
% ELEC 537 Project
% Nate Raymondi and WeiChong Cheng
% MUSE Throughput Simulation for users with bursty traffic
    % NOTE - we have assumed a MIMO AP with rank 4

% INPUTS
    % num_users - number of users connected to the AP
    % prob - VECTOR of probabilities of a packet being available when a
        % user is called upon, should be of length(num_users)
    % trials - number of trials wished to perform, equivalently, this is
        % the number of times a group will send data packets to the AP
        
% OUTPUTS 
    % throughput - returns the average PHY throughput
    % thru - PHY throughput for each transmission attempt, 
        % length(thru_trials) = trials

% Constant overhead values (can be changed, idk what they really are)
t_data = 100e-6;
t_OH = 1e-6;
t_ACK = 5e-9;
t_trig = 100e-9;

% User vector
users = 1:num_users;
% User probability vector - 
%prob = 0:1/num_users:1;
%prob = prob(1:end-1);

% thru_trials = zeros(1,trials);
% for i = 1:trials
    %probability distribution among users
%     prob = rand(1,num_users).*(pr-1)*0.1;
    % random user wins medium
    winner = winner;
    x1 = winner + 1;
    x2 = winner + 2;
    x3 = winner + 3;
    
    % get grouped user indexes, account for user group vector wrapping
    if x1/num_users > 1
        x1 = mod(x1, num_users);
    end
    if x2/num_users > 1
        x2 = mod(x2, num_users);
    end 
    if x3/num_users > 1
        x3 = mod(x3, num_users);
    end
    
    % binary decide if grouped users have packets to send
%     winplus1 = sum(rand >= cumsum([1-prob(x1), 1]));
%     winplus2 = sum(rand >= cumsum([1-prob(x2), 1]));
%     winplus3 = sum(rand >= cumsum([1-prob(x3), 1]));
    winplus0 = 0;
    winplus1 = 0;
    winplus2 = 0;
    winplus3 = 0;
    if(pkts(winner) > 0)
        winplus0 = 1;
        pkts(winner) = pkts(winner) - 1;
    end
    if(pkts(x1) > 0)
        winplus1 = 1;
        pkts(x1) = pkts(x1) - 1;
    end
    if(pkts(x2) > 0)
        winplus2 = 1;
        pkts(x2) = pkts(x2) - 1;
    end
    if(pkts(x3) > 0)
        winplus3 = 1;
        pkts(x3) = pkts(x3) - 1;
    end
    hasData = [winplus0 winplus1 winplus2 winplus3];
    groupMembers = [winner, x1, x2, x3];
    % PHY throughput for this run
%     thru_trials(i) = ((t_data + winplus1*t_data + winplus2*t_data + winplus3*t_data)...
%         /(t_trig + t_data + t_OH + 3*t_ACK));
% end

% Calculate average PHY throughput
% throughput = sum(thru_trials)/trials * 100;
% disp(['groupMember for MUSE are : ',winner,x1,x2,x3]);

end


