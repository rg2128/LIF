ca2 = [];
for kk = 1:16
    for jj =1:2
        for mm =1:4
            for od = 1:2
                for imjk = 1:20
                    a = AllTrace3(imjk,od,mm, jj, kk);
                    a=a{1,1};
                    dt = 0.01;
                    meanfire = [];
                    timebin = [];
                    n_osn = 10;
                    count= 1;
                    i=1;
                    allz = [];
                    while i <length(a)-1
                        start = a(i,1);
                        stop = a(i,1) + dt;
                        if a(i+1,1)>stop
                            meanfire(count) = 1;
                            timebin(count) = a(i+1,1);
                            count = count+1;
                            allz(1) = i+1;
                        else
                            allz = find (a(:,1)>stop);
                            %      allz(1)-i
                            if ~isempty(allz)
                            meanfire(count) = allz(1)-i;
                            timebin(count) = a(allz(1),1);
                            count = count+1;
                            end
                        end
                        if ~isempty(allz)
                            i = allz(1);
                        else
                            break;
                        end
                    end
                    
                    bins = 0.00:0.01:0.5;
                    meanfireNorm = zeros(1,length(bins));
                    for n =1:length(timebin)
                        if timebin(1,n)>0
%                             for p =1:length(bins)
                               first =  find(timebin(1,n) > bins);
%                                 last = find(timebin(1,n) < bins);
                            meanfireNorm(first(end)) = meanfire(1,n);
                            end
                        end
                    
%                     figure;
%                     plot(meanfire)
% 
%                     figure;
%                     plot(meanfire/n_osn)


                    Sim.t_Start=0;
                    Sim.t_End=0.5;
                    % Sim.dt_step=0.01;
                    % dt=Sim.dt_step;            % time step (s)
                    Tseq=Sim.t_Start:dt:Sim.t_End;
                    ca2exp=2;
                    ca2tau=0.15;
                    
                    ca2_filter  = Tseq.*exp(-Tseq/ca2tau);
                    
                    figure;
                    plot(ca2_filter)

                    ca2(imjk,:) = conv((meanfireNorm/n_osn).^ca2exp,ca2_filter);
                    
%                     ca2 = conv((meanfireNorm/n_osn).^ca2exp,ca2_filter);
%                     figure;
%                     plot(ca2)
                end
%                  GoCa2 = PADCAT(ca2{1},);
