
process foo {
    input:
    path(f)
    output:
    tuple val("see"), path("see")
    script:
    """
    wc -l $f > see
    """
}

process bar {
    input:
    tuple val(label), path(f)
    output:
    path("see")
    script:
    """
    wc -l $f > see
    """
}



workflow flow1 {
    take: data
    main:
        foo(data)
        bar(foo.out)
    emit:
        bar.out
}
