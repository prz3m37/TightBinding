import os 
import sys
from mysql.connector import MySQLConnection


class DataBaseConnection(object):

    def __init__(self) -> None:
        self.__host = None
        self.__port = None
        self.__user = None
        self.__password = None
        self.__database = None
        self.__connection = None
        self.__limit_of_reconnections = 0

    def set_connection(self):
        self.__connection = MySQLConnection(
            host=self.__host,
            user=self.__user,
            password=self.__password,
            database=self.__database)
        return

    def close_connection(self):
        self.__connection.close()
        return

    def execute_query(self, query):
        cursor = self.__connection.cursor()
        cursor.execute(query)
        cursor.commit()
        return

    def check_connection(self):
        if self.__connection is None:
            reconnection = 0
            while reconnection <= self.__limit_of_reconnections:
                reconnection += 1
                self.set_connection()
                if self.__connection is not None and self.__connection.is_connected():
                    break
                else:
                    continue
            return
        else:
            return
