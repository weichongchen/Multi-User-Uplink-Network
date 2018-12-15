% function [throughput, thru_trials] = user_selection(num_users, trials, mode, pr)
function [pkts,hasData,groupMembers,ExpTraffic] = user_selection(num_users,winner,pkts,mode,ExpTraffic,succ)
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
t_ACK = 100e-9;
t_trig = 5e-9;
% User vector
users = 1:num_users;

% User probability vector - 
%prob = 0:1/num_users:1;
%prob = prob(1:end-1);

Xs0 = zeros(1,3); %smoothing component
Xt0 = zeros(1,3); %trend component
Xs1 = zeros(1,3);
Xt1 = zeros(1,3);
Xs2 = zeros(1,3);
Xt2 = zeros(1,3);
Xs3 = zeros(1,3);
Xt3 = zeros(1,3);
% thru_trials = zeros(1,trials);
% for i = 1:trials
    %probability distribution among users
%     prob = rand(1,num_users).*(pr-1)*0.1;
%     ExpTraffic = prob;%zeros(1,num_users);
    %for the initial case
    if(sum(ExpTraffic)==0 || sum(ExpTraffic==max(ExpTraffic))==num_users)
        % random user wins medium
        winner = randi([1 num_users],1,1);
        x1 = winner + 1;
        x2 = winner + 2;
        x3 = winner + 3;
    %part-fair--backoff+three highest probability
    else
        winner = winner;
        Traffic_buffer = ExpTraffic;
        Traffic_buffer(winner) = 0;
        [value index] = sort(Traffic_buffer,'descend');
        x1 = index(1);
        x2 = index(2);
        x3 = index(3);
        
    %zero-fair--choose only the four highest probability
%     else
%         [value index] = sort(ExpTraffic,'descend');
%         winner = index(1);
%         x1 = index(2);
%         x2 = index(3);
%         x3 = index(4);
    end
    
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
%     winplus0 = sum(rand >= cumsum([1-prob(winner), 1]));
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
    
    if(mode==1)
        %% simple moving average
        i = succ;
        if(i>1)
            ExpTraffic(winner) = (ExpTraffic(winner)*(i-1) + winplus0)/i;
            ExpTraffic(x1) = (ExpTraffic(x1)*(i-1) + winplus1)/i;
            ExpTraffic(x2) = (ExpTraffic(x2)*(i-1) + winplus2)/i;
            ExpTraffic(x3) = (ExpTraffic(x3)*(i-1) + winplus3)/i;
        else
            ExpTraffic(winner) = (ExpTraffic(winner) + winplus0)/(i+1);
            ExpTraffic(x1) = (ExpTraffic(x1) + winplus1)/(i+1);
            ExpTraffic(x2) = (ExpTraffic(x2) + winplus2)/(i+1);
            ExpTraffic(x3) = (ExpTraffic(x3) + winplus3)/(i+1);
        end
    elseif(mode==2)
        %% weighted moving average
        alpha = 0.5;
        ExpTraffic(winner) = ((1-alpha)*ExpTraffic(winner) + alpha*winplus0);
        ExpTraffic(x1) = ((1-alpha)*ExpTraffic(x1) + alpha*winplus1);
        ExpTraffic(x2) = ((1-alpha)*ExpTraffic(x2) + alpha*winplus2);
        ExpTraffic(x3) = ((1-alpha)*ExpTraffic(x3) + alpha*winplus3);
    else
        %% Holt-Winters
        alpha = 0.5;
        beta = 0.6;
        %winner
        Xf0 = Xs0 + Xt0;
        Xs0(3) = alpha*winplus0 + (1-alpha)*Xf0(2);
        Xt0(3) = beta*(Xs0(2) - Xs0(1)) + (1-beta)*Xt0(1);
        Xf0 = Xs0 + Xt0;
        ExpTraffic(winner) = Xf0(2);
        Xs0(1:2) = Xs0(2:3);
        Xs0(3) = 0;
        Xt0(1:2) = Xt0(2:3);
        Xt0(3) = 0;

        %x1
        Xf1 = Xs1 + Xt1;
        Xs1(3) = alpha*winplus0 + (1-alpha)*Xf1(2);
        Xt1(3) = beta*(Xs1(2) - Xs1(1)) + (1-beta)*Xt1(1);
        Xf1 = Xs1 + Xt1;
        ExpTraffic(x1) = Xf1(2);
        Xs1(1:2) = Xs1(2:3);
        Xs1(3) = 0;
        Xt1(1:2) = Xt1(2:3);
        Xt1(3) = 0;

        %x2
        Xf2 = Xs2 + Xt2;
        Xs2(3) = alpha*winplus0 + (1-alpha)*Xf2(2);
        Xt2(3) = beta*(Xs2(2) - Xs2(1)) + (1-beta)*Xt2(1);
        Xf2 = Xs2 + Xt2;
        ExpTraffic(x2) = Xf2(2);
        Xs2(1:2) = Xs2(2:3);
        Xs2(3) = 0;
        Xt2(1:2) = Xt2(2:3);
        Xt2(3) = 0;

        %x3
        Xf3 = Xs3 + Xt3;
        Xs3(3) = alpha*winplus0 + (1-alpha)*Xf3(2);
        Xt3(3) = beta*(Xs3(2) - Xs3(1)) + (1-beta)*Xt3(1);
        Xf3 = Xs3 + Xt3;
        ExpTraffic(x3) = Xf3(2);
        Xs3(1:2) = Xs3(2:3);
        Xs3(3) = 0;
        Xt3(1:2) = Xt3(2:3);
        Xt3(3) = 0;
    end

    % PHY throughput for this run
%     thru_trials(i) = ((t_data + winplus1*t_data + winplus2*t_data + winplus3*t_data)...
%         /(2*t_trig + t_data + t_OH + 3*t_ACK));
% end

% Calculate average PHY throughput
% throughput = sum(thru_trials)/trials * 100;
% if(mode==1)
%     disp(['groupMember for SMA are : ',winner,x1,x2,x3]);
% elseif(mode==2)
%     disp(['groupMember for WMA are : ',winner,x1,x2,x3]);
% else
%     disp(['groupMember for HW are : ',winner,x1,x2,x3]);
% end

end


