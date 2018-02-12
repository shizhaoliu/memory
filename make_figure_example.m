function make_figure_example(example_results_super,example_results_unsuper,state_idx)
    % Plot supervsied results

    % Plot trajectory of two event
    figure1 = figure(1);
    set(figure1,'Position',[300,10,600,800]);
    subplot(3,2,1)
    imagesc(example_results_super(1).prob')
    colormap(flipud(gray))
    set(gca,'Ydir','Normal')
    xlabel('Time bin')
    ylabel('Position bin')
    title('Example 1')
    subplot(3,2,2)
    imagesc(example_results_super(2).prob')
    colormap(flipud(gray))
    set(gca,'Ydir','Normal')
    xlabel('Time bin')
    title('Example 2')
    suptitle('Supervised results')

    % Plot distribution of shuffled Rw
    % Dash line marks the threshold value of Rw(95%)
    % Solid line marks the real Rw value of this event
    % If solid line is larger than dash line, this event is considered significant
    subplot(3,2,3)
    histogram(abs(example_results_super(1).Rw_null{1}),20,'FaceColor',[1 1 1])
    hold on
    plot([example_results_super(1).thrs_w(1),example_results_super(1).thrs_w(1)],[0,200],'--','color','black','linewidth',2)
    hold on
    plot([abs(example_results_super(1).Rw(1)),abs(example_results_super(1).Rw(1))],[0,200],'color','black','linewidth',2)
    ylabel('Number')
    xlabel('|Rw|')
    [p,~] = SignificanceTest(abs(example_results_super(1).Rw), abs(example_results_super(1).Rw_null{1}));
    if p < 0.001
        title('p < 0.001');
    else
        title(['p=',num2str(p)]);
    end
    subplot(3,2,4)
    histogram(abs(example_results_super(2).Rw_null{1}),20,'FaceColor',[1 1 1])
    hold on
    plot([example_results_super(2).thrs_w(1),example_results_super(2).thrs_w(1)],[0,200],'--','color','black','linewidth',2)
    hold on
    plot([abs(example_results_super(2).Rw(1)),abs(example_results_super(2).Rw(1))],[0,200],'color','black','linewidth',2)
    xlabel('|Rw|')
    [p,~] = SignificanceTest(abs(example_results_super(2).Rw), abs(example_results_super(2).Rw_null{1}));
    if p < 0.001
        title('p < 0.001');
    else
        title(['p=',num2str(p)]);
    end
    % Plot distribution of shuffled Rwd
    % Dash line marks the threshold value of Rwd(95%)
    % Solid line marks the real Rwd value of this event
    % If solid line is larger than dash line, this event is considered significant

    subplot(3,2,5)
    histogram(abs(example_results_super(1).Rwd_null{1}),20,'FaceColor',[1 1 1])
    hold on
    plot([example_results_super(1).thrs_wd(1),example_results_super(1).thrs_wd(1)],[0,200],'--','color','black','linewidth',2)
    hold on
    plot([abs(example_results_super(1).Rwd(1)),abs(example_results_super(1).Rwd(1))],[0,200],'color','black','linewidth',2)
    xlabel('Rwd')
    ylabel('Number')
    [p,~] = SignificanceTest(example_results_super(1).Rwd, example_results_super(1).Rwd_null{1});
    if p < 0.001
        title('p < 0.001');
    else
        title(['p=',num2str(p)]);
    end
    subplot(3,2,6)
    histogram(abs(example_results_super(2).Rwd_null{1}),20,'FaceColor',[1 1 1])
    hold on
    plot([example_results_super(2).thrs_wd(1),example_results_super(2).thrs_wd(1)],[0,200],'--','color','black','linewidth',2)
    hold on
    plot([abs(example_results_super(2).Rwd(1)),abs(example_results_super(2).Rwd(1))],[0,200],'color','black','linewidth',2)
    xlabel('Rwd')
    [p,~] = SignificanceTest(example_results_super(2).Rwd, example_results_super(2).Rwd_null{1});
    if p < 0.001
        title('p < 0.001');
    else
        title(['p=',num2str(p)]);
    end
    %%
    % Plot unsupervsied results

    % Plot trajectory of two event
    figure2 = figure(2);
    set(figure2,'Position',[300,100,600,600]);
    subplot(2,2,1)
    % For visualization of unsupervised result,since the order of state is
    % sorted by occupancy, we need to sort the state based on their according
    % position first. Note that this is only for visualization.
    imagesc(example_results_unsuper(1).prob(:,state_idx)')
    colormap(flipud(gray))
    set(gca,'Ydir','Normal')
    xlabel('Time bin')
    ylabel('State')
    title('Example 1')
    subplot(2,2,2)
    imagesc(example_results_unsuper(2).prob(:,state_idx)')
    colormap(flipud(gray))
    set(gca,'Ydir','Normal')
    xlabel('Time bin')
    title('Example 2')
    suptitle('Unsupervised results')
    % Plot distribution of shuffled Rwd
    % Dash line marks the threshold value of Rwd(95%)
    % Solid line marks the real Rwd value of this event
    % If solid line is larger than dash line, this event is considered significant
    subplot(2,2,3)
    histogram(abs(example_results_unsuper(1).Rwd_null{1}),20,'FaceColor',[1 1 1])
    hold on
    plot([example_results_unsuper(1).thrs_wd(1),example_results_unsuper(1).thrs_wd(1)],[0,200],'--','color','black','linewidth',2)
    hold on
    plot([abs(example_results_unsuper(1).Rwd(1)),abs(example_results_unsuper(1).Rwd(1))],[0,200],'color','black','linewidth',2)
    xlabel('Rwd')
    ylabel('Number')
    [p,~] = SignificanceTest(example_results_unsuper(1).Rwd, example_results_unsuper(1).Rwd_null{1});
    if p < 0.001
        title('p < 0.001');
    else
        title(['p=',num2str(p)]);
    end
    subplot(2,2,4)
    histogram(abs(example_results_unsuper(2).Rwd_null{1}),20,'FaceColor',[1 1 1])
    hold on
    plot([example_results_unsuper(2).thrs_wd(1),example_results_unsuper(2).thrs_wd(1)],[0,200],'--','color','black','linewidth',2)
    hold on
    plot([abs(example_results_unsuper(2).Rwd(1)),abs(example_results_unsuper(2).Rwd(1))],[0,200],'color','black','linewidth',2)
    xlabel('Rwd')
    [p,~] = SignificanceTest(example_results_unsuper(2).Rwd, example_results_unsuper(2).Rwd_null{1});
    if p < 0.001
        title('p < 0.001');
    else
        title(['p=',num2str(p)]);
    end