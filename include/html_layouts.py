class html_layout :
    def __init__( self ) :
        pass
    @staticmethod
    def return_layout_1( ) :
        layout_1 = """
<style>
    body {
        width: 100vw;
        height: 100vh;
        display: flex;
        flex-wrap: wrap;
    }
    .image-container {
        width: 50%;
        display: flex;
        flex-direction: column;
        justify-content: flex-start;
        align-items: center;
        padding-top: 2%;
    }
    .image-title {
        margin-bottom: 10px;
    }
    .table-container {
        width: 100%;
        overflow-x: auto;
    }
    table {
        width: 100%; 
        font-size: 0.8rem;
        border-collapse: collapse;
        margin-top: 20px;
    }
    table, th, td {
        border:1px solid black;
    }
    th, td {
        border: 1px solid black;
        padding: 5px;
        text-align: center;
        white-space: nowrap;
    }
</style>
        """
        return layout_1

    # @staticmethod
    # def return_layout_2( ) :