%% Figure 1 
%Code to plot Figure 1 in King et al., 2024. 

Input_data = importdata('Kingetal2024_IceModel_Inputd18O.txt'); %Input raw Lisiecki & Raymo d18O data

Time = flip(Input_data(:,1)); %Create Time vector from input
d18O = flip(Input_data(:,2)); %Create d18O vector from input

neg_d18O = (d18O)*-1; %Get the negative of the d18O curve
min_d18O = min(neg_d18O); %Create min d18O vector
max_d18O = max(neg_d18O); %Create max d18) vector
Scalar_factor_initial = zeros(length(Time),1); %Create scalar_factor_intitial vector and fill with Time in first column and 0s in 2nd

for i = 1:length(Time) %For the length of the Time vector
    Scalar_factor_initial(i) = (neg_d18O(i)-min_d18O)/(max_d18O-min_d18O); %Normalize d18O curve 
end


%Figure 1a
subplot(2,1,1)
plot(Time,d18O,'color','k','LineWidth',2)
xlim([2900 3650])
ylim([2.8 3.8])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',18)
set(gcf,'color','w')
set(gca,'XDir','reverse')
set(gca,'YDir','reverse')
xlabel('Time (ka)','FontSize',20)
ylabel('?18O (?)','FontSize',20)

%Figure 1b
subplot(2,1,2)
plot(Time,Scalar_factor_initial,'color','k','LineWidth',2)
xlim([2900 3650])
ylim([-0.1 1.2])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',18)
set(gcf,'color','w')
set(gca,'XDir','reverse')
xlabel('Time (ka)','FontSize',20)
ylabel('Normalized Stack)','FontSize',20)