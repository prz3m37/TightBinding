from mysql.connector import MySQLConnection


class DataBaseConnection(object):

    def __init__(self, settings, helpers) -> None:

        self.__connection = None
        self.__helpers = helpers
        self.__host = settings["host"]
        self.__port = settings["port"]
        self.__user = settings["user"]
        self.__table = settings["table"]
        self.__password = settings["password"]
        self.__database = settings["data_base"]

        self.__limit_of_reconnections = 0

    def set_connection(self):
        self.__helpers.save_log(
            '[INFO]: Setting up connection with data base\n')
        self.__connection = MySQLConnection(
            host=self.__host,
            user=self.__user,
            password=self.__password,
            database=self.__database)
        return

    def close_connection(self):
        self.__connection.close()
        self.__helpers.save_log(
            '[INFO]: Closing up connection with data base\n')
        return

    def execute_query(self, query, query_type):
        cursor = self.__connection.cursor()
        cursor.execute(query)
        if query_type == "save":
            cursor.commit()
            return
        else:
            data = cursor.fetchall()
            return data

    def check_connection(self):
        self.__helpers.save_log(
            '[INFO]: Checking if connection with data base exists\n')
        if self.__connection is None:
            self.__helpers.save_log(
                '[INFO]: Restarting connection with data base\n')
            reconnection = 0
            while reconnection <= self.__limit_of_reconnections:
                reconnection += 1
                self.set_connection()
                if self.__connection is not None and self.__connection.is_connected():
                    self.__helpers.save_log(
                        '[ERROR]: Connection with data base has been restarted\n')
                    break
                else:
                    self.__helpers.save_log(
                        '[INFO]: Reconnection number: {}\n'.format(reconnection))
                    continue
            return
        else:
            self.__helpers.save_log('[INFO]: Data base connection exists\n')
            return
