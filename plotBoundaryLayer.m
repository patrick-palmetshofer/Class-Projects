clear
close all

file = dir('output.csv');

table=readtable(file.name,'Delimiter',',');
mat = table.Variables;
figure
hold on
for i = 2:4
    plot(mat(:,1),mat(:,i),'LineWidth',2);
end
grid
set(gca,'FontSize',13)
legend('F','F''','F''''')
xlabel('\eta')
ylabel('F/F''/F''''')