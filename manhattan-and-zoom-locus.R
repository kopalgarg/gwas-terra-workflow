library(shiny)
library(ggplot2)
library(LDlinkR)
library(gaston)
require(BuenColors)
require(biomaRt)
require(qqman)
require(dplyr)
require(data.table)
require(magrittr)
require(gridExtra)
require(ggrastr)

toy_df <- read_tsv('../exampleGWAS.tsv')
ui <- fluidPage(

    titlePanel("GWAS"),
    

    sidebarLayout(
        sidebarPanel(
            fileInput(inputId = 'in_data', label = 'Upload association data', multiple = F),
            selectInput(inputId = 'chr',label = 'CHR',choices = c(1:22), selected = 2),
            textInput(inputId = 'pos', label = 'POS', value = 60718043),
            selectInput(inputId = 'pop', label = 'POP', choices=c('EUR', 'ASN', 
                                                                  'AFR', 'AMR', 
                                                                  'ALL'), selected = 'EUR'),
            tableOutput('significant_table')
        ),

        mainPanel(
           plotOutput("manhattan_plot"),
           plotOutput('zoom_plot'))
    )
)

server <- function(input, output) {

    output$manhattan_plot <- renderPlot({
        data = toy_df
        if (!is.null(input$in_data)){
            data = read_tsv(input$in_data$datapath)
        }
        manhattan_ggplot(data)
    })
    output$significant_table <- renderTable({
        data = toy_df
        if (!is.null(input$in_data)){
            data = read_tsv(input$in_data$datapath)
        }
        significant(data)
    })
    output$zoom_plot <- renderPlot({
        chr = as.numeric(input$chr)
        pos = as.numeric(input$pos)
        pop = input$pop
        data = toy_df
        if (!is.null(input$in_data)){
        data = read_tsv(input$in_data$datapath)
        }
        zoomlocus(data, chr, pos, pop)
    })
    manhattan_ggplot <- function(data){
        
        
        don <- data %>% 
            
            group_by(CHR) %>% 
            summarise(chr_len=as.numeric(max(BP))) %>% 
            
            mutate(tot=cumsum(chr_len)-chr_len) %>%
            select(-chr_len) %>%
            
            left_join(data, ., by=c("CHR"="CHR")) %>%
            
            arrange(CHR, BP) %>%
            mutate( BPcum=BP+tot) %>%
            
            mutate( is_annotate=ifelse(-log10(P)>7, "yes", "no"))
        
        sub =don %>% group_by(CHR) %>% mutate(min =min(P)) %>% filter(min==P) %>% filter(-log10(min)>7.5)
        
        
        axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
        
        p = ggplot(don, aes(x=BPcum, y=-log10(P))) +
            
            ggrastr::geom_point_rast( aes(color=as.factor(CHR)), alpha=0.8, size=0.1) +
            scale_color_manual(values = rep(c("grey", "black"), 22 )) +
            
            scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
            scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
            
            geom_label_repel( data=sub, aes(label=SNP), size=2) +
            geom_hline(yintercept=log10(10e4), linetype="dashed", 
                       color = "blue", size=0.5)+
            
            geom_hline(yintercept=log10(31000000), linetype="dashed", 
                       color = "red", size=0.5)+
            xlab('CHR') +
            
            theme_bw() +
            theme( 
                legend.position="none",
                panel.border = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank()
            )
        p
        
    }
    significant <- function(data){
        
        
        don <- data %>% 
            
            group_by(CHR) %>% 
            summarise(chr_len=as.numeric(max(BP))) %>% 
            
            mutate(tot=cumsum(chr_len)-chr_len) %>%
            select(-chr_len) %>%
            
            left_join(data, ., by=c("CHR"="CHR")) %>%
            
            arrange(CHR, BP) %>%
            mutate( BPcum=BP+tot) %>%
            
            mutate( is_annotate=ifelse(-log10(P)>7, "yes", "no"))
        
        sub =don %>% group_by(CHR) %>% mutate(min =min(P)) %>% filter(min==P) %>% filter(-log10(min)>6) %>%
            subset(., select= c('SNP','CHR','BP','P'))
        return(sub)
        
    }
    zoomlocus <- function(data, sel.chr, sel.pos, pop){
        
        token= '4be7ca3e6929'
        range = 500000
        dat.bmi.sel.region <- data %>% filter(CHR == sel.chr, between(BP, sel.pos - 
                                                                          range, sel.pos + range))
        ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
        out.bm.genes.region <- getBM(
            attributes = c('start_position','end_position','ensembl_gene_id','external_gene_name', 'gene_biotype'), 
            filters = c('chromosome_name','start','end'), 
            values = list(sel.chr, sel.pos - range,sel.pos+ range), 
            mart = ensembl)
        LD_p=LDproxy(paste0('chr',sel.chr, ':', sel.pos), pop = pop, r2d = "r2", token = token)
        LD_p$Coord=as.character(LD_p$Coord)
        LD_p$BP=as.numeric(gsub(LD_p$Coord, pattern = paste0('chr',sel.chr,':'), replacement = ''))
        merged=left_join(dat.bmi.sel.region,LD_p, by='BP') %>% filter(!is.na(pos))
        merged_main = merged %>% filter(BP==as.numeric(sel.pos))
        
        b=ggplot(data = merged, aes(x=BP, y=-log10(P), color=R2)) + geom_point(size=2,alpha=0.7) +
            scale_color_gradientn(colors = jdb_palette("solar_basic"))+ pretty_plot() +
            geom_label(data=merged_main,size=4, aes(label=paste0(CHR,':',BP))) +
            theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank()) +
            theme(legend.position = "none") +
            theme(axis.title.y=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank()) +
            theme(legend.position = "top") 
        
        
        out.bm.genes.region <- getBM(
            attributes = c('start_position','end_position','ensembl_gene_id','external_gene_name', 'gene_biotype'), 
            filters = c('chromosome_name','start','end'), 
            values = list(sel.chr, min(merged$BP), max(merged$BP)), 
            mart = ensembl)
        plot.range=c(min(merged$BP),max(merged$BP))
        p2 <- ggplot(data = out.bm.genes.region) + 
            geom_linerange(aes(x = external_gene_name, ymin = start_position, ymax = end_position, colour = gene_biotype, group = gene_biotype)) +
            coord_flip() + ylab("") +
            ylim(plot.range) + 
            geom_text(aes(x = external_gene_name, y = start_position, label = external_gene_name, colour = gene_biotype), fontface = 2, alpha = I(0.7), hjust = "right", size= 2.5) + 
            expand_limits(y=c(-1, 1)) +
            pretty_plot() +
            theme(axis.title.y=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank()) +
            theme(legend.position = "bottom") 
        
        grid.arrange(b,p2)
    }

}

# Run the application 
shinyApp(ui = ui, server = server)