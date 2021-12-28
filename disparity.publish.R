# libraries
{
library(strucchange)
library(ggplot2)
library(ggpubr)

library(XLConnect)
library(rgl)
library(smooth)
library(forecast)
library(gap)
library(lme4)
library(dendextend)
library(mclust)
library(factoextra)
library(gap)
library(psych)
library(car)
library(reshape2)
}

###################### Define Functions ######################
my_forward = function(data.in, y.col='y', x.col, cv.col, cutoff=0.05) {
	cat('... forward selection ... '); flush.console();
	if(!is.null(cv.col)) cv.col = cv.col[which(!is.na(cv.col))]
	if(is.null(cv.col) || length(cv.col) == 0 || is.na(cv.col) || nchar(cv.col) == 0) cv.col=c()
	cc = c(y.col, x.col, cv.col)
	data.all = data.frame(data.in[, cc])
	colnames(data.all)[1:length(cc)] = cc
	colnames(data.all)[1] = 'y'
	
	col = colnames(data.all)
	len = length(col)
	kept <- c()
	col_checked <- c();
	flag <- T;
	while(flag) {  
		## Find the best single covariate on conditioned response
		target <- 0
		p_target <- 1
#		cat('try... ');
		col.excl = which(col %in% c('y', cv.col, col_checked))
		col.test = 1:len
		col.test = col.test[which(!col.test %in% col.excl)]
		for(i in col.test) {  
			string <- paste("y ~ ", paste(c(cv.col, kept, col[i]), collapse="+"), sep="")
			fit <- lm(data=data.all, formula=as.formula(string), na.action=na.exclude)
			coeff <- summary(fit)$coefficients
			p <- coeff[-c(1:(length(c(cv.col, kept))+1)), 4];
			p_adj <- min(p)  # for categorical variate with multiple categories, find the smallest p
			if(!is.na(p_adj) & p_adj < p_target) { 
				p_target <- p_adj
				target <- i
			}
#			cat(col[i], '..')
		}
		target_name <- col[target];
		col_checked <- c(col_checked, target_name);
#		cat('\ntarget_name:', target_name, '\n')

		if(p_target > cutoff) {
			flag <- F;
		} else {
			if(length(kept) == 0) {
				kept <- target_name
			} else {
				kept_test <- c(kept, target_name)
				string <- paste("y ~ ", paste(c(cv.col, kept_test), collapse="+"), sep="")
				fit <- lm(data=data.all, formula=as.formula(string), na.action=na.exclude)
				coeff <- summary(fit)$coefficients[-c(1:(length(cv.col)+1)), ]
				var_names <- rownames(coeff);
				p <- coeff[ ,4];
				p_adj <- p 
				sig_names <- var_names[which(p_adj <= 0.05)]
				kept <- kept_test[unique(unlist(lapply(kept_test, function(x) grep(x, sig_names))))]
			}
		}
#		cat('kept:', kept, '\n')
		flush.console();
				
		if(length(col_checked) == (len-1-length(cv.col))) {
			flag <- F;
		}
	}           
	cat('selected:', kept, '\n')

	if(!exists('kept') || length(kept) == 0 || is.na(kept)) {
		return(c())
	} else {
		fit <- lm(data=data.all, formula = as.formula(paste("y ~ ", paste(c(cv.col, kept), collapse="+"), sep="")), na.action=na.exclude)
		coeff <- summary(fit)$coefficients
		var_names <- row.names(coeff);
		p <- coeff[, 4];
		p_adj <- p 
		effect <- coeff[, 1];
		result <- data.frame(var=var_names, p=p, eff=effect, stringsAsFactors=F);
		result <- result[-1, ]
		return(result)
	}
}

my_backward = function(data.in, y.col='y', x.col, cv.col, cutoff=0.05) {
	cat('... backward selection ... '); flush.console();
	if(!is.null(cv.col)) cv.col = cv.col[which(!is.na(cv.col))]
	if(is.null(cv.col) || length(cv.col) == 0 || is.na(cv.col) || nchar(cv.col) == 0) cv.col=c()
	cc = c(y.col, x.col, cv.col)
	data.all = data.frame(data.in[, cc])
	colnames(data.all)[1:length(cc)] = cc
	colnames(data.all)[1] = 'y'
	
	col = colnames(data.all)
	len = length(col)
	flag <- T;
	while(flag) {  
		## Find the worst single covariate on conditioned response
		fit = lm(y ~ ., data=data.all)
		coeff = summary(fit)$coefficients[-1, ]
		coeff = coeff[which(!rownames(coeff) %in% cv.col), ]
		if(length(coeff) < 5) coeff = t(data.frame(coeff))
		worst.idx = which.max(coeff[, 4])
		if(length(worst.idx) == 0) {
			flag = F;
		} else if(coeff[worst.idx, 4] > cutoff) {
			worst.col = rownames(coeff)[worst.idx]
			if(nrow(coeff) == 1) {
				flag=F
			} else {
				data.all = data.all[, which(colnames(data.all) != worst.col)]
			}
		} else {
			flag = F
		}
		dim(data.all)		
	}
	kept = colnames(data.all)
	kept = kept[which(!kept %in% c('y', cv.col))]
	cat('selected:', kept, '\n')

	if(!exists('kept') || length(kept) == 0 || is.na(kept)) {
		return(c())
	} else {
		fit <- lm(data=data.all, formula = as.formula(paste("y ~ ", paste(c(cv.col, kept), collapse="+"), sep="")), na.action=na.exclude)
		coeff <- summary(fit)$coefficients
		var_names <- row.names(coeff);
		p <- coeff[, 4];
		p_adj <- p 
		effect <- coeff[, 1];
		result <- data.frame(var=var_names, p=p, eff=effect, stringsAsFactors=F);
		result <- result[-1, ]
		return(result)
	}
}

load('Disparity/dfr.list.RData')

###################### Find Stages ######################
{
	cum = do.call(rbind, lapply(dfr.list, function(x) {
		s = sum(x$case.cnt)
		p = sum(x$pop.max)
		r = c(s, p, s/p*1000, log(s/p*1000, 10))
		return(r)
	}))
	cum = data.frame(cum)
	colnames(cum) = c('case.cnt', 'pop', 'case.per.1k', 'case.norm.log')
	dim(cum)  ## 26  3
	cum$new.case.cnt = c(0, cum$case.cnt[-1] - cum$case.cnt[-nrow(cum)])
	cum$new.case.per.1k = cum$new.case.cnt / cum$pop * 1000
	cum$new.case.norm.log = log(cum$new.case.per.1k, 10)
	cum$date = rownames(cum)
	cum$biweek = 1:nrow(cum)

	par(mfrow=c(2,1))
	start=2; end=nrow(cum); cum.wk = cum[start:end, ]; y = 'new.case.norm.log'; h = 3; ## y = 'new.case.cnt'
	brk = breakpoints(as.formula(paste(y, '~ biweek')), data=cum.wk, h=h)
	plot(x=cum.wk$biweek, cum.wk[, y]); abline(v=brk$breakpoints)
	bb = c(1, brk$breakpoints, nrow(cum.wk))
	bb.p = c()
	for(i in 3:length(bb)) {
		window.start = bb[i-2]
		window.end = bb[i]
		window.break = bb[i-1] - bb[i-2] + 1
		window.df = cum.wk[window.start:window.end, ]
		bb.p = c(bb.p, sctest(as.formula(paste(y, '~ biweek')), data=window.df, type="Chow", point=window.break)$p.value)
	}
	(bb.p = data.frame(brk=brk$breakpoints, p=bb.p))		
}
	

## create period lists 
{
	total.pop = sum(dfr.list[[1]]$pop.max)
	dd.case = data.frame(rbind(c(1,6), c(7, 12), c(13, 19), c(20, 23), c(24, 27), c(28, length(dfr.list)))); colnames(dd.case) = c('start', 'end'); ddx.case=1:length(dfr.list); periods.case=c('Surge-A', 'Decline-A', 'Valley', 'Surge-B', 'Decline-B', 'Surge-C')
	fill.colors = c('gray40', 'gray65', 'gray90', 'gray40', 'gray65', 'gray40')

	dfr.period.list = list()
	wk.list = dfr.list; y = 'case.cnt'; y.norm = 'diff'; dd = dd.case; ddx=ddx.case
	for(i in 1:nrow(dd)) {
		start = dd[i, 1]; if(start > 1) start = start - 1;
		end = dd[i, 2] 
		df.start = wk.list[[start]]
		df.end = wk.list[[end]]
		diff = df.end[, y] - df.start[, y]
		diff = ifelse(diff <= 0, 0.5, diff)
		diff = log(diff/df.start$pop.max * 1000, 10)
		df.start[, y.norm] = diff
		dfr.period.list[[i]] = df.start
	}
	length(dfr.period.list)  ## 6
}

## plot & compare growth rate
{
	wk.list = dfr.list; wk.period.list = dfr.period.list; y = 'case.cnt'; y.norm = 'diff'; dd = dd.case; ddx=ddx.case; periods=periods.case
	
	rate = c()
	for(i in 1:length(wk.period.list)) {
		this = wk.period.list[[i]]
		start = dd[i, 1]; if(start > 1) start = start - 1;
		end = dd[i, 2] 
		rr = this[, y.norm] / (ddx[end] - ddx[start]) / 2
		rate = rbind(rate, cbind(rr, rep(i, length(rr))))
	}
	rate = data.frame(Rate=rate[, 1], Period=factor(periods[rate[, 2]], levels=periods), stringsAsFactors=F)

	compare_means(Rate ~ Period, data=rate, method='t.test')

	fc = c('gray65', 'gray80', 'gray90', 'gray65', 'gray80', 'gray65')
	p = ggboxplot(rate, x='Period', y='Rate', fill='Period') + labs(x='', y='Growth Rate') 
	p = p + scale_fill_manual(values=fc) + theme(legend.position='none') 
	p + stat_compare_means(comparisons=list(c('Surge-A', 'Decline-A'), c('Decline-A', 'Valley'), c('Valley', 'Surge-B'), c('Surge-B', 'Decline-B'), c('Decline-B', 'Surge-C')), method='t.test', label='p.signif', label.y=seq(0.4, 0.65, length.out=length(periods.case)))
}

###################### Feature Selection ######################
flag = 'grow.all';  ## all.flag = 'cnt.few' if total count at the first and the last time point; all.flag = 'cnt.all' if at all time point
{
	features.all = c('pop.ratio.sex', 'pct.age.65.over', 'pct.hisp', 'pct.aind', 'pct.asian', 'pct.black', 'pct.poverty', 'pct.unemploy', 'household.income.med'); ## 
	cv.col = c('pop.density', 'pop.age.median');
	length(features.all)  ## 13
	cutoff = 0.05

	wk.list = dfr.list; wk.period.list = dfr.period.list; y = 'case.cnt'; y.norm = 'diff'; y2 = 'adjustedcasenumbers_log10'; dd = dd.case; ddx=ddx.case
	if(flag == 'cnt.few') {
		y.norm = y2; wk.list = wk.list[c(1, length(wk.list))]; wk.period.list = list(wk.period.list[[1]])
	} else if(flag == 'cnt.all') {
		y.norm = y2;
	}
		
	selected.indv = c()
	for(i in 1:length(wk.list)) {
#			cat('time', i, '...\n'); flush.console()		
		df.wk = wk.list[[i]]
		fw = my_forward(data.in=df.wk, y.col=y.norm, x.col=features.all, cv.col=cv.col, cutoff=cutoff)
		if(!is.null(fw)) {
			bw = my_backward(data.in=df.wk, y.col=y.norm, x.col=fw$var[which(!fw$var %in% cv.col)], cv.col=cv.col, cutoff=cutoff)
			bw$time = i;
			selected.indv = rbind(selected.indv, bw[, c('var', 'time')])
		}
		bw = my_backward(data.in=df.wk, y.col=y.norm, x.col=features.all, cv.col=cv.col, cutoff=cutoff)
		if(!is.null(bw)) {
			fw = my_backward(data.in=df.wk, y.col=y.norm, x.col=bw$var[which(!bw$var %in% cv.col)], cv.col=cv.col, cutoff=cutoff)
			fw$time = i;
			selected.indv = rbind(selected.indv, fw[, c('var', 'time')])
		}
	}
	selected.indv = selected.indv[which(!selected.indv$var %in% cv.col), ]
	selected.indv = unique(selected.indv)
		
	selected.period = c()
	for(i in 1:length(wk.period.list)) {
#			cat('period', i, '...\n'); flush.console()
		df.wk = wk.period.list[[i]]
		fw = my_forward(data.in=df.wk, y.col=y.norm, x.col=features.all, cv.col=cv.col, cutoff=cutoff)
		if(!is.null(fw)) {
			bw = my_backward(data.in=df.wk, y.col=y.norm, x.col=fw$var[which(!fw$var %in% cv.col)], cv.col=cv.col, cutoff=cutoff)
			bw$time = i;
			selected.period = rbind(selected.period, bw[, c('var', 'time')])
		}
		bw = my_backward(data.in=df.wk, y.col=y.norm, x.col=features.all, cv.col=cv.col, cutoff=cutoff)
		if(!is.null(bw)) {
			fw = my_backward(data.in=df.wk, y.col=y.norm, x.col=bw$var[which(!bw$var %in% cv.col)], cv.col=cv.col, cutoff=cutoff)
			fw$time = i;
			selected.period = rbind(selected.period, fw[, c('var', 'time')])
		}
	}
	selected.period = selected.period[which(!selected.period$var %in% cv.col), ]
	selected.period = unique(selected.period)

}	

###################### Test Associations ######################
{
	cv = cv.col;
	wk.list = dfr.list; wk.period.list = dfr.period.list; y = 'case.cnt'; y.norm = 'diff'; y2 = 'adjustedcasenumbers_log10'; dd = dd.case; ddx=ddx.case; periods=c('Surge-A', 'Decline-A', 'Valley', 'Surge-B', 'Decline-B', 'Surge-C')
	if(flag == 'cnt.few') {
		y.norm = y2; wk.list = wk.list[c(1, length(wk.list))]; wk.period.list = list(wk.period.list[[1]])
	} else if(flag == 'cnt.all') {
		y.norm = y2;
	}

	selected.all = c(cv, unique(c(selected.indv$var, selected.period$var)))
		
	## individual time point
	assoc.indv = c()
	for(sel in selected.all) {
		assoc.sel = c()
		for(i in 1:length(wk.list)) {
			df.wk = wk.list[[i]]
			selected = selected.indv[which(selected.indv$time == i), 'var']
			selected = c(selected, cv)
			selected = selected[which(selected != sel)]
			s = paste(y.norm, '~', paste(c(sel, selected), collapse='+'))
			model = lm(as.formula(s), data=df.wk)
			assoc = summary(model)$coeff[2, c(1, 4)]
			assoc = c(assoc, i)
			assoc.sel = rbind(assoc.sel, assoc)
		}
		assoc.sel = data.frame(var=sel, est=assoc.sel[, 1], p=assoc.sel[, 2], time=assoc.sel[, 3], stringsAsFactors=F)
		assoc.sel$Biweek = assoc.sel$time
		assoc.sel$p.adj = p.adjust(assoc.sel$p, method='BH')
		assoc.sel$period = cut(assoc.sel$Biweek, breaks=c(dd$start-1, nrow(assoc.sel)), labels=periods)
		assoc.sel$anno = ifelse(assoc.sel$Biweek %in% dd$start, as.character(assoc.sel$period), '')
		assoc.indv = rbind(assoc.indv, assoc.sel)
	}
		
	# period
	assoc.period = c()
	for(sel in selected.all) {
		assoc.sel = c()
		for(i in 1:length(wk.period.list)) {
			df.wk = wk.period.list[[i]]
			selected = selected.period[which(selected.period$time == i), 'var']
			selected = c(selected, cv)
			selected = selected[which(selected != sel)]
			s = paste(y.norm, '~', paste(c(sel, cv, selected), collapse='+'))
			model = lm(as.formula(s), data=df.wk)
			assoc = summary(model)$coeff[2, c(1, 4)]
			assoc = c(assoc, i)
			assoc.sel = rbind(assoc.sel, assoc)
		}
		assoc.sel = data.frame(var=sel, est=assoc.sel[, 1], p=assoc.sel[, 2], time=assoc.sel[, 3], stringsAsFactors=F)
		assoc.sel$p.adj = p.adjust(assoc.sel$p, method='BH')
		assoc.sel$period = periods[1:length(wk.period.list)]
		assoc.sel = cbind(assoc.sel, dd)
		assoc.period = rbind(assoc.period, assoc.sel)
	}

	assoc.indv.all = assoc.indv
	assoc.period.all = assoc.period
}

###################### Plots in paper ######################
{
library(gridExtra)
library(ggforce)
library(showtext) # font_paths(); font_files()
font_add("Arial", "arial.ttf")

## plot of new counts at all time points
theme_update(plot.title = element_text(hjust = 0.5))
dd = dd.case[, 1]; plot.fn = 'Disparity/plot.new.counts.times.pdf'
vvs = unique(c(cv, selected.indv$var)); vvs.display=vvs;
# vvs = c('pct.white.nhisp', 'pct.hisp', 'pct.black', 'pct.asian', 'pct.aind'); vvs.display = c('White Non-Hispanic', 'Hispanic', 'Black', 'Asian', 'American Indian')
# vvs = c('household.income.med', 'pct.poverty', 'pct.household.middle', 'pct.household.income.200K.more', 'pct.household.disability', 'pct.unemploy', 'pct.household.snap'); vvs.display = c('Median Income', 'Poverty', 'Middle Class', 'High Income', 'Disability', 'Unemployment')
# vvs = c('pop.age.median', 'pct.age.65.over'); vvs.display = c('Median Age', 'Senior')
# vvs = c('pct.highschool', 'pct.college.plus'); vvs.display = c('High School', 'College +')
plot.list = list()
for(i in 1:length(vvs)) {
	vv = vvs[i]
	vv.display = vvs.display[i]
	assoc.indv = assoc.indv.all[which(assoc.indv.all$var == vv), ]
	assoc.indv$anno = gsub('Surge', 'Srg', assoc.indv$anno); 
	assoc.indv$anno = gsub('Decline', 'Dcl', assoc.indv$anno); 
	assoc.indv$FDR = assoc.indv$p.adj
	assoc.period = assoc.period.all[which(assoc.period.all$var == vv), ]
	assoc.period = as.data.frame((lapply(assoc.period,  rep,  each=2)))
	assoc.period$x = assoc.period$start - 0.5
	assoc.period[seq(2, nrow(assoc.period), 2), 'x'] = assoc.period[seq(1, nrow(assoc.period), 2), 'end'] + 0.5
	assoc.period$anno = ''; assoc.period[seq(2, nrow(assoc.period), 2), 'anno'] = '*';
	
	scaling.min=1.2; scaling.max=1.5; scaling.text = scaling.max - 0.2
	max.est = max(c(assoc.indv$est, assoc.period$est), na.rm=T)
	min.est = min(c(assoc.indv$est, assoc.period$est), na.rm=T)
	if(max.est <= 0) {
		max.est = abs(min.est) * 0.1
	} else if (min.est >= 0) {
		min.est = -abs(max.est) * 0.1
	}	
	range.est = max.est - min.est
	if(range.est < 0.03 & max.est < 0.002) {
		scaling.max = 1.5
		scaling.text = scaling.max - 0.3
	} else if(range.est < 0.03 & max.est < 0.02) {
		scaling.max = 1.35
		scaling.text = scaling.max - 0.15
	} else if(range.est < 0.03) {
		scaling.max = 1.5
		scaling.text = scaling.max - 0.4
	} else {
		scaling.max = 1.8
		scaling.text = scaling.max - 0.2
	}
	print(c(range.est, max.est, scaling.max, scaling.text))

	plt = ggplot(assoc.indv, aes(x=Biweek, y=est)) + geom_rect(aes(xmin=Biweek-0.5, xmax=Biweek+0.5, fill=period), ymin=min.est*scaling.min, ymax=max.est*scaling.max, alpha=0.5, color=NA) + theme_bw(base_size=18) + ggtitle(vv.display) + theme(plot.title=element_text(hjust=0.5)) + labs(y='Effect') + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.title=element_text(size=15), legend.text=element_text(size=10)) + scale_x_continuous(breaks=seq(0, nrow(assoc.indv), 4), limits=c(0, nrow(assoc.indv)+1)) + scale_y_continuous(limits=c(min.est*scaling.min, max.est*scaling.max)) + scale_color_gradient2(low='blue', mid='gray', high='red', midpoint=0) + scale_fill_manual(values=fill.colors) + geom_text(data=assoc.indv[which(nchar(assoc.indv$anno) > 1), ], aes(label=anno, x=Biweek-0.3), y=max.est*scaling.text, size=6, hjust ='left') + scale_size_continuous(range=c(2,7))
	plt = plt + geom_point(aes(color=est, size=-log10(FDR))) + geom_point(data=assoc.indv[which(assoc.indv$p.adj < 0.1), ], aes(x=Biweek, y=est, size=-log10(FDR)/50), shape=4, stroke=1.5, color='white') # + geom_point(data=assoc.indv[which(assoc.indv$p.adj >= 0.1 & assoc.indv$p < 0.05), ], aes(x=Biweek, y=est), size=0.5, shape=1, stroke=1.5, color='white')	
	plt = plt + geom_line(data=assoc.period[which(assoc.period$p.adj < 0.1), ], aes(x=x, y=est, group=period, col=est), linetype='solid', size=1.5) + geom_line(data=assoc.period[which(assoc.period$p.adj > 0.1), ], aes(x=x, y=est, group=period, col=est), linetype='dashed', size=1.5) + geom_hline(yintercept=0, linetype=1, size=0.7)
	plt = plt + guides(fill=FALSE)  #  + theme(legend.position='bottom') 
	plot.list[[i]] = plt
}
ml = marrangeGrob(plot.list, nrow=1, ncol=1)
ggsave(plot.fn, width=12, height=6, ml)

## plot of total counts at selected time points
theme_update(plot.title = element_text(hjust = 0.5))
dd = dd.case[, 1]; wk.list = dfr.list[c(1, length(dfr.list))]; y2 = 'adjustedcasenumbers_log10'; 


