using Distributions


#точная функция log(exp(a)+exp(b)) для больших чисел
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



#Априорные распределения
prior1 = Beta(0.05, 1.95) #распределение для P(кольцо | не женат и на улице)
prior2 = Beta(1.95, 0.05) #распределение для P(кольцо | женат и на улице)
println("изначально ожидаем P(кольцо | не женат)=", mean(prior1), " ± ", std(prior1))
println("изначально ожидаем P(кольцо | женат)=", mean(prior2), " ± ", std(prior2))
#ожидаем стандартную ошибку в 15 процентных пунктов у нерепрезентативного опроса https://adamobeng.com/download/FastCheapAccurate.pdf
#и ожидаем 64% женатых среди взрослых на улице https://belstat.gov.by/upload/iblock/471/471b4693ab545e3c40d206338ff4ec9e.pdf
prior3 = Beta(6.5, 3.5)
println("ожидаем P(женат | взрослый на улице)=", mean(prior3), " ± ", std(prior3))

prior4 = Beta(1.23, 0.04) #P(правая рука | обручальное кольцо) основываясь на https://www.stranamam.ru/post/843327/#comments
println("ожидаем P(правая рука | обручальное кольцо)=", mean(prior4), " ± ", std(prior4))
prior5 = Beta(1,1) #P(правая рука | не обручальное кольцо)
println("ожидаем P(правая рука | не обручальное кольцо)=", mean(prior5), " ± ", std(prior5))


#уличные данные
R_ring = 57
R_n = 164
L_ring = 33
L_n = 258
B_ring = 35
B_n = 103
#ringStreet = R_ring + B_ring
#noRingStreet = R_n + B_n - ringStreet

#данные на реддите https://www.reddit.com/r/AskWomen/comments/xw4kp6/married_women_how_often_do_you_wear_your/
ringReddit = 19 + 16
noRingReddit = 9 + 5

#данные stranamam https://www.stranamam.ru/post/843327/#comments
ringWomen = 198 + 6
noRingWomen = 47 + 11

#данные на хабре https://habr.com/ru/articles/782456/
ringHabr = 25 + 1 + 1
noRingHabr = 47

avg1 = 0
avg2 = 0
avg = -Inf
n = -Inf
avg3 = -Inf
n3 = -Inf
#Количество итераций
niter = 1e+6
num=0
#перебираем модели
for i in (1:niter)
    p1 = rand(prior1, 1)[1] #P(кольцо | на улице и не женат)
    p2 = rand(prior2, 1)[1] #P(кольцо | на улице и женат)
    p3 = rand(prior3, 1)[1] #P(женат | взрослый на улице)

    p9 = rand(prior4, 1)[1] #P(правая | обручальное кольцо)
    p10 = rand(prior5, 1)[1] #P(правая | не обручальное кольцо)

    if p1 > 0 && p1 < 1 && p2 > 0 && p2 < 1 && p3 > 0 && p3 < 1 && p9>0 && p9<1 && p10>0 && p10<1

        #нерепрезентативные оценки P(кольцо | женат) с ste=0.15
        BdivideA = (1.0 / p2) - 1.0 #beta/alpha for beta distribution
        AplusB = p2 * p2 * BdivideA / (0.15 * 0.15) #beta+alpha for beta distribution
        beta = AplusB * BdivideA / (BdivideA + 1)
        alpha = AplusB - beta
        prior6 = Beta(alpha, beta)
        p6 = rand(prior6, 1)[1] #P(кольцо | женатый и на реддите)
        p7 = rand(prior6, 1)[1] #P(кольцо | женатый и на хабре)
        p8 = rand(prior6, 1)[1] #P(кольцо | женатый и на женском форуме)

        p4 = (1 - p3) * (1 - p1) / ((1 - p3) * (1 - p1) + p3 * (1 - p2)) #P(не женат | взрослый на улице без кольца)
        p44 = p3 * p2 / (p3*p2 + (1-p3)*p1) #P(женат | взрослый на улице c кольцом)

        p51 = p3 * p2 * p9 + (1 - p3) * p1 *p10  #P(кольцо на правой | на улице)
        p52 = p3 * p2 *(1-p9) + (1 - p3) * p1*(1-p10) #P(кольцо на левой  | на улице)
        global avg1 += p4
        global avg2 +=p44

        #правдоподобность P(кольцо| на улице) исходя из наблюдений на улице
        loglikelihood = log(p51) * R_ring + log(1 - p51) * (R_n - R_ring) 
        loglikelihood += log(p52) * L_ring + log(1 - p52) * (L_n - L_ring) 
        loglikelihood += log(p51+p52) * R_ring + log(1 - p51 - p52) * (R_n - R_ring) 
        
        
        loglikelihood += log(p6) * ringReddit + log(1 - p6) * noRingReddit #правдоподобность P(кольцо | женат) исходя из наблюдений на реддите
        loglikelihood += log(p7) * ringHabr + log(1 - p7) * noRingHabr #правдоподобность P(кольцо | женат) исходя из наблюдений на хабре
        loglikelihood += log(p8) * ringWomen + log(1 - p8) * noRingWomen #правдоподобность P(кольцо | женат) исходя из наблюдений на stranamam


        global n = logsum(n, loglikelihood)
        global avg = logsum(avg, loglikelihood + log(p4))
        global avg3 = logsum(avg3, loglikelihood + log(p44))

        global num+=1
    end
end

println("изначально ожидали P(не женат | взрослый на улице без кольца) = ", avg1 / num)
println("исходя из данных P(не женат | взрослый на улице без кольца) = ", exp(avg - n))

println("изначально ожидали P(женат | взрослый на улице с кольцом) = ", avg2 / num)
println("исходя из данных P(женат | взрослый на улице с кольцом) = ", exp(avg3 - n))
