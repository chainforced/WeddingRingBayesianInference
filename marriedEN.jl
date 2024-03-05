using Distributions


#precise function log(exp(a)+exp(b)) for big numbers
function logsum(a, b)
    if isnan(a)
        a = -Inf
    end
    if isnan(b)
        b = -Inf
    end
    if isinf(a) && isinf(b)
        return -Inf
    end
    mx = maximum([a, b])
    mn = a + b - mx
    #(a/b+1)*b = a+b
    return log(exp(mn - mx) + 1) + mx
end



#prior distributions
prior1 = Beta(0.05, 1.95) #prior distr for P(ring | not married & on the street)
prior2 = Beta(1.95, 0.05) #prior distr for P(ring | married & on the street)
println("prior expectation: P(ring | not married)=", mean(prior1), " ± ", std(prior1))
println("prior expectation: P(ring | married)=", mean(prior2), " ± ", std(prior2))
#nonrepresentative surveys have STE=0.15 https://adamobeng.com/download/FastCheapAccurate.pdf
#There's 64% married people among adults in Belarus  https://belstat.gov.by/upload/iblock/471/471b4693ab545e3c40d206338ff4ec9e.pdf
prior3 = Beta(6.5, 3.5)
println("expectation: P(married | взрослый на улице)=", mean(prior3), " ± ", std(prior3))

prior4 = Beta(1.23, 0.04) #P(right hand | wedding ring) based on https://www.stranamam.ru/post/843327/#comments
println("expectation: P(right hand | wedding ring)=", mean(prior4), " ± ", std(prior4))
prior5 = Beta(1,1) #P(right hand | not wedding ring)
println("expectation: P(right hand | not wedding ring)=", mean(prior5), " ± ", std(prior5))


#Street Data
R_ring = 57
R_n = 164
L_ring = 33
L_n = 258
B_ring = 35
B_n = 103

#Reddit Data https://www.reddit.com/r/AskWomen/comments/xw4kp6/married_women_how_often_do_you_wear_your/
ringReddit = 19 + 16
noRingReddit = 9 + 5

#stranamam Data https://www.stranamam.ru/post/843327/#comments
ringWomen = 198 + 6
noRingWomen = 47 + 11

#habrahabr data  https://habr.com/ru/articles/782456/
ringHabr = 25 + 1 + 1
noRingHabr = 47

avg1 = 0
avg2 = 0
std1 = 0
std2 = 0
avg = -Inf
n = -Inf
avg3 = -Inf
n3 = -Inf
std0=-Inf
std3=-Inf
#number of iterations for montecarlo simulation
niter = 1e+6
num=0
for i in (1:niter)
    p1 = rand(prior1, 1)[1] #P(ring | not married & on the street)
    p2 = rand(prior2, 1)[1] #P(кольцо | married & on the street)
    p3 = rand(prior3, 1)[1] #P(женат | adult & on the street)

    p9 = rand(prior4, 1)[1] #P(right hand | wedding ring)
    p10 = rand(prior5, 1)[1] #P(right hand | not a wedding ring)

    if p1 > 0 && p1 < 1 && p2 > 0 && p2 < 1 && p3 > 0 && p3 < 1 && p9>0 && p9<1 && p10>0 && p10<1

        #nonrepresentative estimations of P(ring | married) with STE=0.15
        BdivideA = (1.0 / p2) - 1.0 #beta/alpha for beta distribution
        AplusB = p2 * p2 * BdivideA / (0.15 * 0.15) #beta+alpha for beta distribution
        beta = AplusB * BdivideA / (BdivideA + 1)
        alpha = AplusB - beta
        prior6 = Beta(alpha, beta)
        p6 = rand(prior6, 1)[1] #P(ring | married && on reddit)
        p7 = rand(prior6, 1)[1] #P(ring | married && on habr)
        p8 = rand(prior6, 1)[1] #P(ring | married && on women forum)

        p4 = (1 - p3) * (1 - p1) / ((1 - p3) * (1 - p1) + p3 * (1 - p2)) #P(not married | adult on the street without a wedding ring)
        p44 = p3 * p2 / (p3*p2 + (1-p3)*p1) #P(married | adult on the street with a wedding ring)

        p51 = p3 * p2 * p9 + (1 - p3) * p1 *p10  #P(ring on right hand | on the street)
        p52 = p3 * p2 *(1-p9) + (1 - p3) * p1*(1-p10) #P(ring on left hand  | on the street)
        global avg1 += p4
        global avg2 +=p44
        global std1 += p4*p4
        global std2 +=p44*p44

        #loglikelihood of P(ring| on the streed) based on street data
        loglikelihood = log(p51) * R_ring + log(1 - p51) * (R_n - R_ring) 
        loglikelihood += log(p52) * L_ring + log(1 - p52) * (L_n - L_ring) 
        loglikelihood += log(p51+p52) * R_ring + log(1 - p51 - p52) * (R_n - R_ring) 
        
        
        loglikelihood += log(p6) * ringReddit + log(1 - p6) * noRingReddit #loglikelihood of P(ring | married) based on reddit data
        loglikelihood += log(p7) * ringHabr + log(1 - p7) * noRingHabr #loglikelihood of P(ring | married) based on habrahabr data
        loglikelihood += log(p8) * ringWomen + log(1 - p8) * noRingWomen #loglikelihood of P(ring | married) based on women forum data


        global n = logsum(n, loglikelihood)
        global avg = logsum(avg, loglikelihood + log(p4))
        global avg3 = logsum(avg3, loglikelihood + log(p44))
        global std0 = logsum(std0, loglikelihood + log(p4)*2)
        global std3 = logsum(std3, loglikelihood + log(p44)*2)

        global num+=1
    end
end

println("prior expectation: P(не женат | adult on the street without a wedding ring) = ", avg1 / num," ± ", (std1/num - avg1*avg1/(num*num))^0.5 )
println("posterior: P(не женат | adult on the street without a wedding ring) = ", exp(avg - n)," ± ", (exp(std0 - n) - exp(2*(avg - n)))^0.5)

println("prior expectation: P(married | adult on the street with a wedding ring) = ", avg2 / num," ± ",  (std2/num - avg2*avg2/(num*num))^0.5 )
println("posterior: P(married | adult on the street with  awedding ring) = ", exp(avg3 - n)," ± ", (exp(std3 - n) - exp(2*(avg3 - n)))^0.5)
