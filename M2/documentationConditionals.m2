isEven = i -> i % 2 == 0

oddOrEven = i -> (
    if isEven i then "even"
    else "odd")
