% Author: Ankit Kaushik
% Docmented created: 07.06.2013
% Last Modication: 05.10.2014
% Purpose: The script implements reinforcement learning on the number of
% channels N, where it is allowed to sense only M slot, the aim is to
% maximize the long term capacity over a given number of slots
% Following stratergies are implemented to select the channels to be sensed
% in a given slot:
% 1) Random search
% 2) based on utilization probability
% 3) based on belief vector which is sun of the uti probability and a penalty
% penalty = E[To] = 1/alpha is given when a channel is selected and was found busy

clc;
clear all;
close all;
simulation = 0;                                                       % enable this to run the simulation otherwise the res 
if simulation
    %% Declaring variables
    Kilosec = 1e3;                                                    % Kilo second = 1000
    N = [4];                                                         % Total number of channels to be sensed
    alpha_vec = [];                                                   % The probability of starting of disturbance   
    beta_vec = [];                                                    % The probability of ending of disturbance
    u_vec = [];                                                       % Channel utilization vector, for given beta and alpha calcultes the spectrum occupancy probability                                     
    F = [];                                                           % Channel state for i_th time slot
    F_scaled = [];                                                    % Channels at ith time slot are scaled to represent unit in mini slots
    T_slot = 4;                                                      % Number of mini slots comprising of one complete time slot,                                                           
    M = 2;
    M_vec = 0:1:M;                                                    % Number of mini slots used for sensing 
    sim_length = 1e5;                                                 % Number of time slots to simulate, one time unit is one TDMA frame (1 msec)
                                                                      % time slot unit will be scaled to T_slot so as to represent new time unit in mini slot 
                                                                      % defined only for sensing, where M will be sensing interval  
    sim_length_secs = sim_length/Kilosec; 

    sensing_sequence_ref_vec = [];                                    % Sensing sequence of the channels corresponding to the ref case (random)
    sensing_sequence_uti_vec = [];                                    % Sensing sequence of the channels corresponding to the utilization case
    sensing_sequence_uti_pen_vec = [];                                % Sensing sequence of the channels corresponding to the channel utilization with penalty case

    flag_opp_found_ref = 0;                                           % Flag to indicate that first oppurtunity found by ref (random) will be used as home channel   
    flag_opp_found_uti = 0;                                           % Flag to indicate that first oppurtunity found by uti seq will be used as home channel
    flag_opp_found_uti_pen = 0;                                       % Flag to indicate that oppurtunity found by uti_pen seq seq will be used as home channel  

    Time_scale = [1:1:sim_length_secs];                               % Time axis denoting the time in seconds 1000 slots = 1 sec (1 sec)
    oppurtunities_ref = zeros(1,length(M_vec));                       % Ratio of oppurtunities found per sec / total number of oppurtunities exiting per sec for ref seq
    oppurtunities_uti = zeros(1,length(M_vec));                       % Ratio of oppurtunities found per sec / total number of oppurtunities exiting per sec for uti seq
    oppurtunities_uti_pen = zeros(1,length(M_vec));                   % Ratio of oppurtunities found per sec / total number of oppurtunities exiting per sec for utilization with penalty seq
    oppurtunities_idle = zeros(1,length(M_vec));                      % Ratio of oppurtunities found per sec / total number of oppurtunities exiting per sec for idle seq

    oppurtunities_ref_th = zeros(1,length(M_vec));                    % Ratio of oppurtunities found per sec / total number of oppurtunities exiting per sec for ref seq determined theoretically
    oppurtunities_uti_th = zeros(1,length(M_vec));                    % Ratio of oppurtunities found per sec / total number of oppurtunities exiting per sec for uti seq determined theoretically
    oppurtunities_uti_pen_th = zeros(1,length(M_vec));                % Ratio of oppurtunities found per sec / total number of oppurtunities exiting per sec for utilization with penalty seq determined theoretically
    oppurtunities_idle_th = zeros(1,length(M_vec));                   % Ratio of oppurtunities found per sec / total number of oppurtunities exiting per sec for idle seq determined theoretically
    
    cap_ref = zeros(1,length(M_vec));                                 % Capacity in bits per sec for ref seq
    cap_uti = zeros(1,length(M_vec));                                 % Capacity in bits per sec for uti seq
    cap_uti_pen = zeros(1,length(M_vec));                             % Capacity in bits per sec for utilization with penalty seq
    cap_idle = zeros(1,length(M_vec));                                % Capacity in bits per sec under idle policy  
    
    total_opputunities = zeros(1,length(M_vec));                      % Total number of oppurtunities popped up by channel in a given time slot
    total_opputunities_th = zeros(1,length(M_vec));                   % Total number of oppurtunities popped up by channel in a given time slot determined theoretically

    
    penalty_vec = [];                                                 % When a channel is found busy will be given a penalty in oder to learn from the past                                     
    belief_vec = [];                                                  % Results of penalty will be stored inside the belief vector   

    penalty_wght_vec = [];                                            % When a channel is found busy will be given a penalty in oder to learn from the past + assigned weight                                    
    belief_wght_vec = [];                                             % Results of penalty and weights will be stored inside the belief vector 

    chan_real = 1e0;

    for i = 1:chan_real
        
        disp(strcat('chan real = ',num2str(i)));         
        randseq_alpha = randperm(N,N);                               % Choice of alpha corresponding to each channel is idle
        randseq_beta = randperm(N,N);                                % Choice of beta corresponding to each channel is also distributed randomly
        alpha_vec = 0.5*rand(1,N);% 0.3*randseq_alpha/N                  % [0.8*1/N(i):0.8*1/N(i):0.8]; 
        beta_vec = 0.3*rand(1,N); %0.3*randseq_alpha/N                   % 0.2*ones(1,N(i));            

        channel_utilization_vec = beta_vec./(alpha_vec + beta_vec);
        
        %Analysis                
        mean_penalty = 1./alpha_vec;
        var_penalty = (1-alpha_vec)./(alpha_vec.^2);

        % Thoretical calculation of oppurtunities
        total_opputunities_th = sum(channel_utilization_vec);
        
        % Set of Primary channels
        F = zeros(1,N);

        % Random sequence does not follow any sensing sequence

        % Utilization case: Is always sorted according to the utilization sequence 
        [B sensing_sequence_uti_vec] = sort(channel_utilization_vec,'ascend');   %[1:1:N(i)];

        % Utilization + penalty: Initial sorting seq is ascending order of the
        % utilization and then updated according to the sensing sequence
        [B sensing_sequence_uti_pen_vec] = sort(channel_utilization_vec,'ascend');    

        % Penalty when a channel is found busy during sensing
        penalty_vec = zeros(1, N);
        belief_vec = zeros(1, N);

        penalty_wght_vec = zeros(1, N);
        belief_wght_vec = zeros(1, N);

        F_scaled = zeros(N,T_slot); 

        % Number of simulated slots
        for i=1:sim_length
            % Update sensing sequence ref
            %sensing_sequence_ref_vec = circshift(sensing_sequence_ref_vec, [1 (length(sensing_sequence_ref_vec) + 1) - home_channel_ref]);
            sensing_sequence_ref_vec = randperm(N,N);

            %% Learn from history : Utilization + penalty case
            % Update the penalties for the sequence  
            for l=1:N
                if penalty_vec(l) ~= 0
                    penalty_vec(l) = penalty_vec(l) - 1;                   
                end
            end

            % Update sensing sequence uti_pen
            % find the index of the new home channel
            % Assign the lowest weight i.e. highest priorty
            weights = ones(1,length(sensing_sequence_uti_pen_vec));

            % Updating the belief vector
            for l=1:N
                channel = sensing_sequence_uti_pen_vec(l);
                belief_vec(channel) = penalty_vec(channel) + weights(channel) * channel_utilization_vec(channel); 
            end 

            % sorting the sensing sequence depending on the
            % priorities and lowest weight assigned to home channel
            % will be given highest priority as the sequence is
            % sorted based on the home channel wieght and         
            [B sensing_sequence_uti_pen_vec] = sort(belief_vec,'ascend'); 
            F = channels_state_for_next_time_slot(N, beta_vec, alpha_vec, F);

            F_scaled = zeros(N,T_slot);
            % Scaling of the time slots into mini slots (sensing intervals)
            F_scaled = F'*ones(1,T_slot);       

            %% Finding idle channel in the sensing interval i.e. with M_vec
            for j=1:length(M_vec)
                % Increment total oppurtunities for the time slot
                total_opputunities(j) = total_opputunities(j) + length(find (F == 1));
                
                % Sensing based on idle criterion 
                oppurtunities_idle(j) = oppurtunities_idle(j) + min(M_vec(j), length(find (F == 1)));

                for k=1:M_vec(j)

                    % Sensing based on ref criterion 
                    if F_scaled(sensing_sequence_ref_vec(k),k) == 1
                        % Idle channel found
                        % Inc the oppurtunity count for the sec
                        oppurtunities_ref(j) = oppurtunities_ref(j) + 1; 
                        % determine sensing overhead ref
                    end

                    % Sensing based on uti criterion 
                    if F_scaled(sensing_sequence_uti_vec(k),k) == 1
                        % Idle channel found
                        % Inc the oppurtunity count for the sec
                        oppurtunities_uti(j) = oppurtunities_uti(j) + 1; 
                        % determine sensing overhead ref

                        % If the oppurtunity found is the first than home
                        % channel will be updated and for the rest of the
                        % opprutunuites on the rest of the channels will only be
                        % detected but will be not set to home channel
                        if flag_opp_found_uti == 0
                            % Update home channel
                            home_channel_uti = sensing_sequence_uti_vec(j);
                            flag_opp_found_uti = 1;
                        end      
                    end
                    % Sensing based on uti_pen criterion
                    if F_scaled(sensing_sequence_uti_pen_vec(k),k) == 1

                        % Idle channel found
                        % Inc the oppurtunity count for the sec
                        oppurtunities_uti_pen(j) = oppurtunities_uti_pen(j) + 1; 
                        % determine sensing overhead uti_pen

                        % If the oppurtunity found is the first than home
                        % channel will be updated and for the rest of the
                        % opprutunuites on the rest of the channels will only be
                        % detected but will be not set to home channel
                        if flag_opp_found_uti_pen == 0
                            % Update home channel
                            home_channel_uti_pen = sensing_sequence_uti_pen_vec(j);
                            flag_opp_found_uti_pen = 1;
                        end               

                    else
                        %% Learn, when the oppurtunity is not detected than it recieves a penalty
                        if (alpha_vec(sensing_sequence_uti_pen_vec(k)) < 0.5) 
                            penalty_vec(sensing_sequence_uti_pen_vec(k)) = round(1/(alpha_vec(sensing_sequence_uti_pen_vec(k))));  
                        end
                    end
                end
            end        
        end
    end
    for l = 1:length(M_vec)        
        cap_ref(l) = oppurtunities_ref(l)/total_opputunities(l);
        cap_uti(l) = oppurtunities_uti(l)/total_opputunities(l);
        cap_uti_pen(l) = oppurtunities_uti_pen(l)/total_opputunities(l);
        cap_idle(l) = oppurtunities_idle(l)/total_opputunities(l);
    end
%    save('results_RLearning.mat');
%    exit;
end
load('results_RLearning.mat');
y = [cap_ref(2), cap_uti(2), cap_uti_pen(2), cap_idle(2)] * 100
bar(y,0.4);
set(gca,'XTick',1:4,'XTickLabel', {'Random', 'Static Learning', 'R. Learning', 'Optimum'})
pbaspect([1 0.62 1]);
Fontsize = 9;
axis tight;
ylabel('Exploited Opportunities $\%$', 'FontSize',Fontsize);
set(gca,'FontSize',Fontsize);
laprint(1, 'Learning', 'options', 'factory', 'width', 8, 'scalefonts',...
    'on', 'factor',0.5, 'keepfontprops', 'on');