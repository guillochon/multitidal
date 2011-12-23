subroutine read_table_dims(filename, rows, cols)
    implicit none
    integer u
    integer, intent(out) :: rows, cols
    character(50), intent(in) :: filename
    parameter (u=21)
    open(u, file=filename, status='old')
    read(u,*) rows, cols
    close(u)
    return
end

subroutine read_table(filename, table, rows, cols)
    implicit none
    integer i, j, u
    integer, intent(in) :: rows, cols
    real, intent(inout), dimension(rows, cols) :: table
    real dummy
    character(50), intent(in) :: filename
    parameter (u=20)

    open(u, file=filename, status='old')

    read(u, *) dummy, dummy
    read(u, *) ((table(i, j),j=1,cols),i=1,rows)

    close(u)
    return
end
